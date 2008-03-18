module Output_module

  use Logging_module  
  
  implicit none

  private

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petsclog.h"

  PetscInt, parameter :: TECPLOT_INTEGER = 0
  PetscInt, parameter :: TECPLOT_REAL = 1

  PetscInt, parameter :: TECPLOT_FILE = 0
  PetscInt, parameter ::  HDF5_FILE = 1

  PetscInt, parameter :: LIQUID_PHASE = 1
  PetscInt, parameter :: GAS_PHASE = 2

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  public :: Output, OutputTecplot, OutputHDF5, OutputVectorTecplot, &
            OutputBreakthrough, OutputGetVarFromArray
  
contains

! ************************************************************************** !
!
! Output: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine Output(realization,plot_flag)

  use Realization_module
  use Option_module
  
  implicit none
  
  type(realization_type) :: realization
  logical :: plot_flag

  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr)

  ! check for plot request from active directory
  if (.not.plot_flag) then

    if (realization%option%use_touch_options) then
      word = 'plot'
      if (OptionCheckTouch(word)) then
        realization%output_option%plot_name = 'plot'
        plot_flag = .true.
      endif
    endif

  endif

  if (plot_flag) then
  
    if (realization%output_option%print_hdf5) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_hdf5, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call OutputHDF5(realization)
      call PetscLogEventEnd(logging%event_output_hdf5, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr) 
      if (realization%option%myrank == 0) &
        print *, '      Seconds to write to HDF5 file: ', (tend-tstart)
    endif
   
    if (realization%output_option%print_tecplot) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_tecplot, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call OutputTecplot(realization)
      call PetscLogEventEnd(logging%event_output_tecplot, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr) 
      if (realization%option%myrank == 0) &
        print *, '      Seconds to write to Tecplot file(s): ', (tend-tstart)
    endif
  
    realization%output_option%plot_number = realization%output_option%plot_number + 1
  
  endif
  
  call OutputBreakthrough(realization)

  if (realization%option%compute_statistics) then
    call ComputeFlowCellVelocityStats(realization)
    call ComputeFlowFluxVelocityStats(realization)
    call ComputeFlowMassBalance(realization)
  endif
  
  plot_flag = .false.

  call PetscLogStagePop(ierr)
  
end subroutine Output

! ************************************************************************** !
!
! Output_Breakthrough: Main driver for all breakthrough output subroutines
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine OutputBreakthrough(realization)
                           ! for some flakey reason, current Intel 10.1 reports
                           ! error if 'only' statement not used.
  use Realization_module, only : realization_type 
  use Option_Module
  
  implicit none
  
  type(realization_type) :: realization

#if 0  
  if (realization%output_option%print_hdf5) then
    call OutputBreakthroughHDF5(realization)
  endif
#endif
 
  if (realization%output_option%print_tecplot) then
    call OutputBreakthroughTecplot(realization)
  endif

end subroutine OutputBreakthrough


! ************************************************************************** !
!
! OutputTecplot: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputTecplot(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  use Mphase_module
  use Richards_module
  use Richards_Lite_module
  
  use Reactive_Transport_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec
  Vec :: natural_vec
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.tec'
    output_option%plot_name = ''
  else
    if (output_option%plot_number < 10) then
      write(filename,'("pflotran00",i1,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(filename,'("pflotran0",i2,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(filename,'("pflotran",i3,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(filename,'("pflotran",i4,".tec")') output_option%plot_number  
    endif
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit

    ! initial portion of header
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]"'

    ! write flow variables
    string2 = ''
    select case(option%iflowmode)
      case (MPH_MODE)
        string2 = MphaseGetTecplotHeader(realization)
      case(RICHARDS_MODE)
        string2 = RichardsGetTecplotHeader(realization)
      case(RICHARDS_LITE_MODE)
        string2 = RichardsLiteGetTecplotHeader(realization)
    end select
    string = trim(string) // trim(string2)
    
    ! write transport variables
    if (option%ntrandof > 0) then
      string2 = RTGetTecplotHeader(realization)
      string = trim(string) // trim(string2)
    endif

    ! write material ids
    if (associated(patch%imat)) then
      string = trim(string) // ',"Material_ID"'
    endif

    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    if (realization%discretization%itype == STRUCTURED_GRID) then
      ! count vars in string
      quote_count = 0
      comma_count = 0
      do i=1,len_trim(string)
        ! 34 = '"'
        if (iachar(string(i:i)) == 34) then
          quote_count = quote_count + 1
        ! 44 = ','
        else if (iachar(string(i:i)) == 44 .and. mod(quote_count,2) == 0) then
          comma_count = comma_count + 1
        endif
      enddo
      ! there are comma_count+1 variables
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4)') &
                   option%time/output_option%tconv,grid%structured_grid%nx+1, &
                   grid%structured_grid%ny+1,grid%structured_grid%nz+1
      write(string2,'(i5)') comma_count+1
      string = trim(string) // ', VARLOCATION=([4-' // &
               trim(adjustl(string2)) // ']=CELLCENTERED)'
    else
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4)') &
                   option%time/output_option%tconv,grid%structured_grid%nx, &
                   grid%structured_grid%ny,grid%structured_grid%nz 
    endif
    string = trim(string) // ', DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)
  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)  

  ! write out coordinates
  if (realization%discretization%itype == STRUCTURED_GRID) then
    call WriteTecplotStructuredGrid(IUNIT3,realization)
  else
    call GetCoordinates(grid,global_vec,X_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Y_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Z_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
  endif

  select case(option%iflowmode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)

      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! pressure
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
     ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid mole fractions
          do i=1,option%nspec
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select
  
      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nspec
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select 
#if 0      
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call OutputGetVarFromArray(realization,global_vec,VOLUME_FRACTION,ZERO_INTEGER)
        call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      endif
#endif    
      ! phase
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
      end select
      
    case default
  
  end select
  
  if (option%ntrandof > 0) then
    do i=1,option%ntrandof
      call OutputGetVarFromArray(realization,global_vec,TOTAL_CONCENTRATION,i)
      call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
    enddo
  endif  
  
  ! material id
  if (associated(patch%imat)) then
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
  endif

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  close(IUNIT3)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplot(realization)
  endif
  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       X_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           Z_DIRECTION)
      end select
    endif
  endif
      
end subroutine OutputTecplot

! ************************************************************************** !
!
! OutputVelocitiesTecplot: Print velocities to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesTecplot(realization)
 
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
  character(len=MAXWORDLENGTH) :: filename
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
  
  ! open file
  if (output_option%plot_number < 10) then
    write(filename,'("pflotran_vel00",i1,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(filename,'("pflotran_vel0",i2,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(filename,'("pflotran_vel",i3,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(filename,'("pflotran_vel",i4,".tec")') output_option%plot_number  
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",' // &
             '"vlx [m/y]",' // &
             '"vly [m/y]",' // &
             '"vlz [m/y]"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/y]",' // &
               '"vgy [m/y]",' // &
               '"vgz [m/y]"'
    endif
        if (associated(patch%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    if (realization%discretization%itype == STRUCTURED_GRID) then
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4,'','')') &
                   option%time/output_option%tconv, &
                   grid%structured_grid%nx+1,grid%structured_grid%ny+1,grid%structured_grid%nz+1
      string = trim(string) // ' DATAPACKING=BLOCK'
      if (option%nphase > 1) then
        string = trim(string) // ', VARLOCATION=([4-10]=CELLCENTERED)'
      else
        string = trim(string) // ', VARLOCATION=([4-7]=CELLCENTERED)'
      endif
    else
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4,'','')') &
                   option%time/output_option%tconv, &
                   grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
      string = trim(string) // ' DATAPACKING=BLOCK'
    endif
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)    

  ! write out coorindates
  if (realization%discretization%itype == STRUCTURED_GRID) then
    call WriteTecplotStructuredGrid(IUNIT3,realization)
  else
    call GetCoordinates(grid,global_vec,X_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Y_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Z_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
  endif
  
  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  if (option%nphase > 1) then
    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
  endif

  ! material id
  if (associated(patch%imat)) then
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
  endif
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  close(IUNIT3)
  
end subroutine OutputVelocitiesTecplot

! ************************************************************************** !
!
! OutputFluxVelocitiesTecplot: Print intercellular fluxes to Tecplot file in 
!                              BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputFluxVelocitiesTecplot(realization,iphase, &
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
  
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: local_size, global_size
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: i, j, k
  PetscInt :: local_id, ghosted_id
  PetscInt :: adjusted_size
  PetscInt :: count, iconn
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: array(:)
  PetscInt, allocatable :: indices(:)
  Vec :: global_vec, global_vec2
  PetscReal :: sum, average, max, min , std_dev
  PetscInt :: max_loc, min_loc

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
    
  nullify(array)

  call PetscLogEventBegin(logging%event_output_write_flux_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                          
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  ! open file
  filename = 'pflotran_'
  
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
  
  if (output_option%plot_number < 10) then
    write(string,'("00",i1,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(string,'("0",i2,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(string,'(i3,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(string,'(i4,".tec")') output_option%plot_number  
  endif
  
  filename = trim(filename) // trim(string)
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity flux output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
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
        string = trim(string) // ' vlx [m/y]"'
      case(Y_DIRECTION)
        string = trim(string) // ' vly [m/y]"'
      case(Z_DIRECTION)
        string = trim(string) // ' vlz [m/y]"'
    end select 
    
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    select case(direction)
      case(X_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx-1,grid%structured_grid%ny,grid%structured_grid%nz 
      case(Y_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny-1,grid%structured_grid%nz 
      case(Z_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz-1
    end select 
  
  
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  
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
      if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlyz
        nx_local = grid%structured_grid%nlx-1
      endif
    case(Y_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nx*grid%structured_grid%nz
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlxz
        ny_local = grid%structured_grid%nly-1
      endif
    case(Z_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nxy
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
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
        indices(count) = i+grid%structured_grid%nxs+(j-1+grid%structured_grid%nys)*nx_global+ &
                         (k-1+grid%structured_grid%nzs)*nx_global*ny_global
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
        array(count) = grid%x(grid%nL2G(local_id))
        if (direction == X_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dx(local_id)
      enddo
    enddo
  enddo
  ! warning: adjusted size will be changed in ConvertArrayToNatural
  ! thus, you cannot pass in local_size, since it is needed later
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
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
        array(count) = grid%y(grid%nL2G(local_id))
        if (direction == Y_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dy(local_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
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
        array(count) = grid%z(grid%nL2G(local_id))
        if (direction == Z_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dz(local_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option) 
  call VecZeroEntries(global_vec,ierr)
  call VecGetArrayF90(global_vec,vec_ptr,ierr)
  
  ! place interior velocities in a vector
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
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
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  close(IUNIT3)

  call PetscLogEventEnd(logging%event_output_write_flux_tecplot, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  
end subroutine OutputFluxVelocitiesTecplot

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
  use Patch_module
  
  implicit none

  character(len=MAXWORDLENGTH) :: filename
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

  call PetscLogEventBegin(logging%event_output_vec_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  discretization => realization%discretization
  
  ! open file
  if (option%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
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
    if (associated(patch%imat)) then
      string = trim(string) // ',"Material_ID"'
    endif
    write(fid,'(a)') trim(string)
  
    ! write zone header
    if (realization%discretization%itype == STRUCTURED_GRID) then
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4,'','')') &
                   option%time/realization%output_option%tconv, &
                   grid%structured_grid%nx+1,grid%structured_grid%ny+1,grid%structured_grid%nz+1
      string = trim(string) // ' DATAPACKING=BLOCK'
                                                ! 4=dataset name, 5=material_id
      string = trim(string) // ', VARLOCATION=([4-5]=CELLCENTERED)'
    else
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                   &'', K='',i4,'','')') &
                   option%time/realization%output_option%tconv, &
                   grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
      string = trim(string) // ' DATAPACKING=BLOCK'
    endif
    write(fid,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)  
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)    

  ! write out coorindates
  call GetCoordinates(grid,global_vec,X_COORDINATE)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Y_COORDINATE)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Z_COORDINATE)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call DiscretizationGlobalToNatural(discretization,vector,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  if (associated(patch%imat)) then
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_INTEGER)
  endif
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  close(fid)

  call PetscLogEventEnd(logging%event_output_vec_tecplot, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                            
end subroutine OutputVectorTecplot

! ************************************************************************** !
!
! WriteTecplotDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec(fid,realization,vec,datatype)

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
  
end subroutine WriteTecplotDataSetFromVec

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

  call PetscLogEventBegin(logging%event_output_str_grid_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                              
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  
  if (option%myrank == 0) then
    ! x-dir
    count = 0
    do k=1,nz+1
      do j=1,ny+1
        temp_real = grid%structured_grid%origin(X_DIRECTION)
        write(fid,'(es11.4,x)',advance='no') temp_real
        do i=1,nx
          temp_real = temp_real + grid%structured_grid%dx0(i)
          write(fid,'(es11.4,x)',advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo
    ! y-dir
    count = 0
    do k=1,nz+1
      temp_real = grid%structured_grid%origin(Y_DIRECTION)
      do i=1,nx+1
        write(fid,'(es11.4,x)',advance='no') temp_real
        count = count + 1
        if (mod(count,10) == 0) then
          write(fid,'(a)') ""
          count = 0
        endif
      enddo
      do j=1,ny
        temp_real = temp_real + grid%structured_grid%dy0(j)
        do i=1,nx+1
          write(fid,'(es11.4,x)',advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo
    ! z-dir
    count = 0
    temp_real = grid%structured_grid%origin(Z_DIRECTION)
    do i=1,(nx+1)*(ny+1)
      write(fid,'(es11.4,x)',advance='no') temp_real
      count = count + 1
      if (mod(count,10) == 0) then
        write(fid,'(a)') ""
        count = 0
      endif
    enddo
    do k=1,nz
      temp_real = temp_real + grid%structured_grid%dz0(k)
      do j=1,ny+1
        do i=1,nx+1
          write(fid,'(es11.4,x)',advance='no') temp_real
          count = count + 1
          if (mod(count,10) == 0) then
            write(fid,'(a)') ""
            count = 0
          endif
        enddo
      enddo
    enddo

  endif

  call PetscLogEventEnd(logging%event_output_str_grid_tecplot, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                            
end subroutine WriteTecplotStructuredGrid

! ************************************************************************** !
!
! WriteTecplotDataSet: Writes data from an array within a block
!                      of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSet(fid,realization,array,datatype,size_flag)

  use Realization_module
  use Grid_module
  use Option_module
  use Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscReal :: array(:)
  PetscInt, save :: max_local_size_saved = -1
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i
  PetscMPIInt :: iproc, recv_size
  PetscInt :: max_local_size, local_size
  PetscInt :: istart, iend, num_in_array
  PetscInt :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call PetscLogEventBegin(logging%event_output_write_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
  
  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                       PETSC_COMM_WORLD,ierr)
    local_size = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                         PETSC_COMM_WORLD,ierr)
      max_local_size_saved = max_local_size
      if (option%myrank == 0) print *, 'max_local_size_saved: ', max_local_size
    endif
    max_local_size = max_local_size_saved
    local_size = grid%nlmax
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size+10))
    allocate(integer_data_recv(max_local_size))
    do i=1,local_size
      integer_data(i) = int(array(i))
    enddo
  else
    allocate(real_data(max_local_size+10))
    allocate(real_data_recv(max_local_size))
    do i=1,local_size
      real_data(i) = array(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (option%myrank == 0) then
    if (datatype == TECPLOT_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size) exit
        iend = istart+9
        write(fid,'(10(i3,x))') integer_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      integer_data(1:local_size-iend) = integer_data(iend+1:local_size)
      num_in_array = local_size-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size) exit
        iend = istart+9
        write(fid,'(10(es11.4,x))') real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size-iend) = real_data(iend+1:local_size)
      num_in_array = local_size-iend
    endif
    do iproc=1,option%commsize-1
      call MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data_recv,recv_size,MPI_INTEGER,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          integer_data(num_in_array+1:num_in_array+recv_size) = &
                                             integer_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,'(10(i3,x))') integer_data(istart:iend)
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size,MPI_DOUBLE_PRECISION,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          real_data(num_in_array+1:num_in_array+recv_size) = &
                                             real_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,'(10(es11.4,x))') real_data(istart:iend)
        enddo
        if (iend > 0) then
          real_data(1:num_in_array-iend) = real_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      endif
    enddo
    ! Print the remaining values, if they exist
    if (datatype == 0) then
      if (num_in_array > 0) &
        write(fid,'(10(i3,x))') integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(fid,'(10(es11.4,x))') real_data(1:num_in_array)
    endif
  else
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size,MPI_INTEGER,ZERO_INTEGER,local_size, &
                    PETSC_COMM_WORLD,ierr)
    else
      call MPI_Send(real_data,local_size,MPI_DOUBLE_PRECISION,ZERO_INTEGER,local_size, &
                    PETSC_COMM_WORLD,ierr)
    endif
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

  call PetscLogEventEnd(logging%event_output_write_tecplot, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

end subroutine WriteTecplotDataSet

! ************************************************************************** !
!
! OutputBreakthroughTecplot: Print to breakthrough data to TECPLOT file
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine OutputBreakthroughTecplot(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Patch_module
  use Breakthrough_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: fid, icell
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  type(breakthrough_type), pointer :: breakthrough
  logical, save :: first = .true.

  call PetscLogEventBegin(logging%event_output_breakthrough, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  breakthrough => patch%breakthrough%first
  
  if (associated(breakthrough)) then

    if (option%myrank < 10) then
      write(filename,'("breakthrough_",i1,".tec")') option%myrank  
    else if (option%myrank < 100) then
      write(filename,'("breakthrough_",i2,".tec")') option%myrank  
    else if (option%myrank < 1000) then
      write(filename,'("breakthrough_",i3,".tec")') option%myrank  
    else if (option%myrank < 10000) then
      write(filename,'("breakthrough_",i4,".tec")') option%myrank  
    else if (option%myrank < 100000) then
      write(filename,'("breakthrough_",i5,".tec")') option%myrank  
    endif
  
    ! open file
    fid = 86
    if (first) then
      first = .false.
      open(unit=fid,file=filename,action="write",status="replace")
      ! write header
      ! write title
      write(fid,'(a)',advance="no") '"Time[' // trim(output_option%tunit) // ']"'
      do 
        if (.not.associated(breakthrough)) exit
        do icell=1,breakthrough%region%num_cells
          call WriteBreakthroughHeaderForCell(fid,realization, &
                                              breakthrough%region,icell)
        enddo
        breakthrough => breakthrough%next
      enddo
      write(fid,'(a)',advance="yes") ""
      breakthrough => patch%breakthrough%first
    else
      open(unit=fid,file=filename,action="write",status="old", &
           position="append")
    endif
  
    write(fid,'(1es12.4)',advance="no") option%time/output_option%tconv
    do 
      if (.not.associated(breakthrough)) exit
      do icell=1,breakthrough%region%num_cells
        call WriteBreakthroughDataForCell(fid,realization, &
                                          breakthrough%region%cell_ids(icell))
      enddo
      breakthrough => breakthrough%next
    enddo
    write(fid,'(a)',advance="yes") ""
    close(fid)

  endif

  call PetscLogEventEnd(logging%event_output_breakthrough, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      
end subroutine OutputBreakthroughTecplot

! ************************************************************************** !
!
! WriteBreakthroughHeaderForCell: Print a header for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteBreakthroughHeaderForCell(fid,realization,region,icell)

  use Realization_module
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module
  use Region_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscInt :: icell
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: cell_id_string
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  
  patch => realization%patch
  option => realization%option
  field => realization%field
  grid => patch%grid
  
  write(cell_id_string,*) grid%nL2A(region%cell_ids(icell))
  cell_id_string = trim(region%name) // ' ' //adjustl(cell_id_string)

  select case(option%iflowmode)
    case (MPH_MODE)
      string = ',"X [m] '// trim(cell_id_string) // '",' // &
               '"Y [m] '// trim(cell_id_string) // '",' // &
               '"Z [m] '// trim(cell_id_string) // '",' // &
               '"T [C] '// trim(cell_id_string) // '",' // &
               '"P [Pa] '// trim(cell_id_string) // '",' // &
               '"sl '// trim(cell_id_string) // '",' // &
               '"sg '// trim(cell_id_string) // '",' // &
               '"Ul '// trim(cell_id_string) // '",' // &
               '"Ug '// trim(cell_id_string) // '",'
      do i=1,option%nspec
        write(string2,'(''"Xl('',i2,'') '// trim(cell_id_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,option%nspec
        write(string2,'(''"Xg('',i2,'') '// trim(cell_id_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
#if 0      
      if (option%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction '// trim(cell_id_string) // '"'
      endif
#endif      
      string = trim(string) // ',"Phase '// trim(cell_id_string) // '"'
    case(RICHARDS_MODE,RICHARDS_LITE_MODE)
      if (option%iflowmode == RICHARDS_MODE) then
        string = ',"X [m] '// trim(cell_id_string) // '",' // &
                 '"Y [m] '// trim(cell_id_string) // '",' // &
                 '"Z [m] '// trim(cell_id_string) // '",' // &
                 '"T [C] '// trim(cell_id_string) // '",' // &
                 '"P [Pa] '// trim(cell_id_string) // '",' // &
                 '"sl '// trim(cell_id_string) // '",' // &
                 '"Ul '// trim(cell_id_string) // '"' 
      else
        string = ',"X [m] '// trim(cell_id_string) // '",' // &
                 '"Y [m] '// trim(cell_id_string) // '",' // &
                 '"Z [m] '// trim(cell_id_string) // '",' // &
                 '"P [Pa] '// trim(cell_id_string) // '",' // &
                 '"sl '// trim(cell_id_string) // '"'
      endif
      if (option%iflowmode == RICHARDS_MODE) then
        do i=1,option%nspec
          write(string2,'('',"Xl('',i2,'') '// trim(cell_id_string) // '"'')') i
          string = trim(string) // trim(string2)
        enddo
      endif
#if 0      
      if (option%rk > 0.d0) then
        string = trim(string) // ',"Volume Fraction '// trim(cell_id_string) // '"'
      endif
#endif      
    case default
      string = ',"X [m]",' // &
               '"Y [m]",' // &
               '"Z [m]",' // &
               '"T [C]",' // &
               '"P [Pa]",' // &
               '"sl",' // &
               '"C [mol/L]"'
#if 0               
      if (option%rk > 0.d0) then
        string = trim(string) // ',"Volume Fraction"'
      endif
#endif      
  end select
  write(fid,'(a)',advance="no") trim(string)

end subroutine WriteBreakthroughHeaderForCell

! ************************************************************************** !
!
! WriteBreakthroughDataForCell: Print data for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteBreakthroughDataForCell(fid,realization,local_id)

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
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

100 format(es13.6)
101 format(i)
110 format(',',es13.6)
111 format(',',i)

  ghosted_id = grid%nL2G(local_id)
  ! write out coorindates
  write(fid,110,advance="no") grid%x(ghosted_id)
  write(fid,110,advance="no") grid%y(ghosted_id)
  write(fid,110,advance="no") grid%z(ghosted_id)
  
  select case(option%iflowmode)
    case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)

      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,TEMPERATURE,ZERO_INTEGER,local_id)
      end select

      ! pressure
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,PRESSURE,ZERO_INTEGER,local_id)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,LIQUID_SATURATION,ZERO_INTEGER,local_id)
      end select

      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas saturation
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,GAS_SATURATION,ZERO_INTEGER,local_id)
      end select
    
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid energy
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,LIQUID_ENERGY,ZERO_INTEGER,local_id)
      end select
    
      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas energy
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCell(realization,GAS_ENERGY,ZERO_INTEGER,local_id)
      end select

      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid mole fractions
          do i=1,option%nspec
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCell(realization,LIQUID_MOLE_FRACTION,i-1,local_id)
          enddo
      end select
  
      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nspec
            write(fid,110,advance="no") &
              OutputGetVarFromArrayAtCell(realization,GAS_MOLE_FRACTION,i-1,local_id)
          enddo
      end select 
#if 0      
      ! Volume Fraction
      if (option%rk > 0.d0) then
        write(fid,110,advance="no") &
          OutputGetVarFromArrayAtCell(realization,VOLUME_FRACTION,ZERO_INTEGER,local_id)
      endif
#endif    
      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,RICHARDS_MODE)
          write(fid,111,advance="no") &
            int(OutputGetVarFromArrayAtCell(realization,PHASE,ZERO_INTEGER,local_id))
      end select
      
    case default
  
  end select

end subroutine WriteBreakthroughDataForCell

! ************************************************************************** !
!
! OutputHDF5: Print to HDF5 file in Tecplot compatible format
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
  
#ifndef USE_HDF5
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'write to an HDF5 format.'
    print *
  endif
  stop

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
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(output_option_type), pointer :: output_option
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXWORDLENGTH) :: filename = "pflow.h5"
  character(len=MAXSTRINGLENGTH) :: string
  logical, save :: first = .true.
  PetscReal, pointer :: array(:)
  PetscInt :: i
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdf5_err)
#endif
  if (.not.first) call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id, &
                                      hdf5_err,prop_id)
  if (hdf5_err < 0 .or. first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                     prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    if (option%myrank == 0) print *, '--> creating hdf5 output file: ', filename
  else
    if (option%myrank == 0) print *, '--> appending to hdf5 output file: ', &
                                   filename
  endif
  
  if (first) then

    ! create a group for the coordinates data set
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

!GEH - Structured Grid Dependence - Begin
    ! write out coordinates in x, y, and z directions
    string = "X [m]"
    allocate(array(grid%structured_grid%nx))
    do i=1,grid%structured_grid%nx
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dx0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dx0(i-1)+grid%structured_grid%dx0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nx,array,grp_id)
    deallocate(array)
  
    string = "Y [m]"
    allocate(array(grid%structured_grid%ny))
    do i=1,grid%structured_grid%ny
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dy0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dy0(i-1)+grid%structured_grid%dy0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%ny,array,grp_id)
    deallocate(array)
  
    string = "Z [m]"
    allocate(array(grid%structured_grid%nz))
    do i=1,grid%structured_grid%nz
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dz0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dz0(i-1)+grid%structured_grid%dz0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nz,array,grp_id)
    deallocate(array)
!GEH - Structured Grid Dependence - End

    call h5gclose_f(grp_id,hdf5_err)
    
  endif

  ! create a group for the data set
  write(string,'(''Time:'',es12.4,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)   

  select case(option%iflowmode)
  
    case(MPH_MODE,RICHARDS_MODE, &
         RICHARDS_LITE_MODE)

      ! temperature
      select case(option%iflowmode)
        case (MPH_MODE,RICHARDS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          string = "Temperature"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! pressure
      select case(option%iflowmode)
        case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          string = "Pressure"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          string = "Liquid Saturation"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)  
      end select

      ! gas saturation
      select case(option%iflowmode)
        case (MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          string = "Gas Saturation"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select
      
      ! liquid energy
      select case(option%iflowmode)
        case (MPH_MODE,RICHARDS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          string = "Liquid Energy"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas energy
      select case(option%iflowmode)
        case (MPH_MODE)    
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          string = "Gas Energy"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
    
      ! liquid mole fractions
      select case(option%iflowmode)
        case (MPH_MODE,RICHARDS_MODE)
          do i=1,option%nspec
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i-1)
            write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
      
      ! gas mole fractions
      select case(option%iflowmode)
        case (MPH_MODE)      
          do i=1,option%nspec
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i-1)
            write(string,'(''Gas Mole Fraction('',i4,'')'')') i
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
#if 0    
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call OutputGetVarFromArray(realization,global_vec,VOLUME_FRACTION,ZERO_INTEGER)
        string = "Volume Fraction"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      endif
#endif    
      ! phase
      select case(option%iflowmode)
        case (MPH_MODE)
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          string = "Phase"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
      end select
  
    case default

  end select
  
  if (output_option%print_hdf5_velocities) then

    ! velocities
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
      string = "Gas X-Velocity"
      call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
  
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
      string = "Gas Y-Velocity"
      call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
  
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
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
  
  ! call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  first = .false.
#endif
end subroutine OutputHDF5

#ifdef USE_HDF5
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

  logical, save :: first = .true.
  logical, save :: trick_flux_vel_x = .false.
  logical, save :: trick_flux_vel_y = .false.
  logical, save :: trick_flux_vel_z = .false.

  Vec :: global_vec

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
    
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
  if (first) then
    nx_local = grid%structured_grid%nlx
    ny_local = grid%structured_grid%nly
    nz_local = grid%structured_grid%nlz
    if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (i == 0) trick_flux_vel_x = .true.
    if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (j == 0) trick_flux_vel_y = .true.
    if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (k == 0) trick_flux_vel_z = .true.
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
      if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = .true.
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = .true.
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
        nz_local = grid%structured_grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = .true.
  end select  
  allocate(array(nx_local*ny_local*nz_local))


  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option) 
  call VecZeroEntries(global_vec,ierr)
  call VecGetArrayF90(global_vec,vec_ptr,ierr)
  
  ! place interior velocities in a vector
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
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
        array(count) = vec_ptr(iphase+(direction-1)*option%nphase+ &
                               3*option%nphase*(local_id-1))
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
  
  call VecDestroy(global_vec,ierr)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * output_option%tconv

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                        nx_global,ny_global,nz_global, &
                        nx_local,ny_local,nz_local, &
                        grid%structured_grid%nxs,grid%structured_grid%nys,grid%structured_grid%nzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = .false.
  first = .false.

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
  
  call PetscLogEventBegin(logging%event_output_coordinates_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                            
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
  if (option%myrank == 0) then
    call PetscLogEventBegin(logging%event_h5dwrite_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)     
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_output_coordinates_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 

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
! ConvertArrayToNatural: Converts an array  to natural ordering
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine ConvertArrayToNatural(indices,array, &
                                 local_size,global_size)

  implicit none
  
  PetscInt :: local_size, global_size
  PetscInt :: indices(:)
  PetscReal, pointer :: array(:)
  
  Vec :: natural_vec
  PetscInt, allocatable :: indices_zero_based(:)
  PetscReal, pointer :: vec_ptr(:)
  
  call VecCreate(PETSC_COMM_WORLD,natural_vec,ierr)
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
! OutputGetVarFromArrayAtCell: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function OutputGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)

  use Realization_module
  use Option_module

  use Richards_module, only : RichardsGetVarFromArrayAtCell
  use Richards_Lite_module, only : RichardsLiteGetVarFromArrayAtCell

  implicit none
  
  PetscReal :: OutputGetVarFromArrayAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: local_id

  select case(realization%option%iflowmode)
    case(RICHARDS_MODE)
      OutputGetVarFromArrayAtCell = RichardsGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
    case(RICHARDS_LITE_MODE)
      OutputGetVarFromArrayAtCell = RichardsLiteGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
  end select

end function OutputGetVarFromArrayAtCell

! ************************************************************************** !
!
! OutputGetVarFromArray: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  use Richards_module, only : RichardsGetVarFromArray
  use Richards_Lite_module, only : RichardsLiteGetVarFromArray
  use Reactive_Transport_module, only : RTGetVarFromArray

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  type(option_type), pointer :: option

  call PetscLogEventBegin(logging%event_output_get_var_from_array, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                        
  option => realization%option

  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION,LIQUID_ENERGY, &
         GAS_ENERGY,LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,VOLUME_FRACTION, &
         PHASE,LIQUID_DENSITY,GAS_DENSITY)
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          call RichardsGetVarFromArray(realization,vec,ivar,isubvar)
          return
        case(RICHARDS_LITE_MODE)
          call RichardsLiteGetVarFromArray(realization,vec,ivar,isubvar)
          return
      end select
    case(FREE_ION_CONCENTRATION,TOTAL_CONCENTRATION)
      call RTGetVarFromArray(realization,vec,ivar,isubvar)
    case(MATERIAL_ID)
      if (option%nflowdof > 0) then
        select case(option%iflowmode)
          case(RICHARDS_MODE)
            call RichardsGetVarFromArray(realization,vec,ivar,isubvar)
            return
          case(RICHARDS_LITE_MODE)
            call RichardsLiteGetVarFromArray(realization,vec,ivar,isubvar)
            return
        end select
      else
        call RTGetVarFromArray(realization,vec,ivar,isubvar)
      endif
    case default
      call printErrMsg(realization%option,'IVAR not found in OutputGetVarFromArray')
  end select

  call PetscLogEventEnd(logging%event_output_get_var_from_array, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  
end subroutine OutputGetVarFromArray

! ************************************************************************** !
!
! GetCellCenteredVelocities: Computers the cell-centered velocity component 
!                            as an averages of cell face velocities
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetCellCenteredVelocities(realization,vec,iphase,direction)

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
  PetscInt :: iconn
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: velocity, area
  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, allocatable :: sum_area(:)
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set

  call PetscLogEventBegin(logging%event_output_get_cell_vel, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                            
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
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      area = cur_connection_set%area(iconn)* &
             dabs(cur_connection_set%dist(direction,iconn))
      velocity = patch%internal_velocities(iphase,iconn)* &
                 area
      if (local_id_up > 0) then
        vec_ptr(local_id_up) = vec_ptr(local_id_up) + velocity
        sum_area(local_id_up) = sum_area(local_id_up) + area
      endif
      if (local_id_dn > 0) then
        vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + velocity
        sum_area(local_id_dn) = sum_area(local_id_dn) + area
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%flow_boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      area = cur_connection_set%area(iconn)* &
             dabs(cur_connection_set%dist(direction,iconn))
      vec_ptr(local_id) = vec_ptr(local_id)+ &
                          patch%boundary_velocities(1,iconn)* &
                          area
      sum_area(local_id) = sum_area(local_id) + area
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

  call PetscLogEventEnd(logging%event_output_get_cell_vel, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 

end subroutine GetCellCenteredVelocities

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
  PetscInt :: iconn, i
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
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set

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
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
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
  boundary_condition => patch%flow_boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      vec_ptr(local_id) = vec_ptr(local_id)+ &
                          patch%boundary_velocities(1,iconn)* &
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

  call VecSum(mass_vec,sum,ierr)
  average = sum/real(grid%nmax)
  call VecSet(global_vec,average,ierr)
  call VecMax(mass_vec,max_loc,max,ierr)
  call VecMin(mass_vec,min_loc,min,ierr)
  call VecAYPX(global_vec,-1.d0,mass_vec,ierr)
  call VecNorm(global_vec,NORM_2,std_dev,ierr)
  string = 'Mass Balance'
  if (option%myrank == 0) then
    write(*,'(/,a,/, &
                 &"Average:",1es12.4,/, &
                 &"Max:    ",1es12.4,"  Location:",i11,/, &
                 &"Min:    ",1es12.4,"  Location:",i11,/, &
                 &"Std Dev:",1es12.4,/)') trim(string), &
                                          average,max,max_loc+1, &
                                          min,min_loc+1,std_dev
    write(IUNIT2,'(/,a,/, &
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
  PetscInt :: iconn, i, direction, iphase
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
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set

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
      connection_list => grid%internal_connection_list
      cur_connection_set => connection_list%first
      do 
        if (.not.associated(cur_connection_set)) exit
        do iconn = 1, cur_connection_set%num_connections
          ghosted_id_up = cur_connection_set%id_up(iconn)
          ghosted_id_dn = cur_connection_set%id_dn(iconn)
          local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
          local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
          ! velocities are stored as the downwind face of the upwind cell
          flux = patch%internal_velocities(iphase,iconn)* &
                   cur_connection_set%area(iconn)* &
                   dabs(cur_connection_set%dist(direction,iconn))
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
      boundary_condition => patch%flow_boundary_conditions%first
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection
        do iconn = 1, cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          vec_ptr(local_id) = vec_ptr(local_id)+ &
                              patch%boundary_velocities(iphase,iconn)* &
                              cur_connection_set%area(iconn)* &
                              dabs(cur_connection_set%dist(direction,iconn))
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

      if (option%myrank == 0) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(IUNIT2,'(/,a,/, &
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
  
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: iphase
  PetscInt :: direction
  PetscInt :: local_id, ghosted_id
  PetscInt :: iconn
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, global_vec2
  PetscReal :: sum, average, max, min , std_dev
  PetscInt :: max_loc, min_loc

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
    
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
      connection_list => grid%internal_connection_list
      cur_connection_set => connection_list%first
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
      if (option%myrank == 0) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(IUNIT2,'(/,a,/, &
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

end module Output_module

    
       
