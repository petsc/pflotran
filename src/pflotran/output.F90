module Output_module

  use Logging_module  
  
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petsclog.h"

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
  PetscTruth :: observation_first
  PetscTruth :: hdf5_first
  PetscTruth :: mass_balance_first

  public :: OutputInit, Output, OutputVectorTecplot, &
            OutputObservation, OutputGetVarFromArray

contains

! ************************************************************************** !
!
! OutputInit: Initializes variables
! author: Glenn Hammond
! date: 01/22/09
!
! ************************************************************************** !
subroutine OutputInit(realization)

  use Realization_module

  implicit none
  
  type(realization_type) :: realization
  
  ! set size to -1 in order to re-initialize parallel communication blocks
  max_local_size_saved = -1
  observation_first = PETSC_TRUE
  hdf5_first = PETSC_TRUE
  mass_balance_first = PETSC_TRUE

end subroutine OutputInit

! ************************************************************************** !
!
! Output: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine Output(realization,plot_flag,transient_plot_flag)

  use Realization_module, only : realization_type
  use Option_module, only : OptionCheckTouch, option_type, &
                            output_option_type, printMsg
  
  implicit none
  
  type(realization_type) :: realization
  PetscTruth :: plot_flag
  PetscTruth :: transient_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  type(option_type), pointer :: option

  option => realization%option

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
      call PetscLogEventBegin(logging%event_output_hdf5, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call OutputHDF5(realization)
      call PetscLogEventEnd(logging%event_output_hdf5, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr)
      write(option%io_buffer,'(f6.2," Seconds to write HDF5 file.")') tend-tstart
      call printMsg(option)
    endif
   
    if (realization%output_option%print_tecplot) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_tecplot, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
      if (realization%output_option%tecplot_format == TECPLOT_BLOCK_FORMAT) then
        call OutputTecplotBlock(realization)
      else
        call OutputTecplotPoint(realization)
      endif
      call PetscLogEventEnd(logging%event_output_tecplot, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f6.2," Seconds to write to Tecplot file(s)")') &
            tend-tstart
      call printMsg(option)        
    endif

    if (realization%output_option%print_vtk) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_vtk, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
      call OutputVTK(realization)

      call PetscLogEventEnd(logging%event_output_vtk, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f6.2," Seconds to write to VTK file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
      
    if (realization%output_option%print_mad) then
      call PetscGetTime(tstart,ierr) 
      call PetscLogEventBegin(logging%event_output_mad, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
      call OutputMAD(realization)

      call PetscLogEventEnd(logging%event_output_mad, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      call PetscGetTime(tend,ierr) 
      write(option%io_buffer,'(f6.2," Seconds to write to MAD HDF5 file(s)")') &
            tend-tstart
      call printMsg(option) 
    endif
      
    if (option%compute_statistics) then
      call ComputeFlowCellVelocityStats(realization)
      call ComputeFlowFluxVelocityStats(realization)
!      call OutputMassBalance(realization)
!      call ComputeFlowMassBalance(realization)
    endif
  
    realization%output_option%plot_number = realization%output_option%plot_number + 1

  endif
  
  if (transient_plot_flag) then
    if (option%compute_mass_balance_new) then
      call OutputMassBalanceNew(realization)
    endif
    call OutputObservation(realization)
  endif
  
  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  realization%output_option%plot_name = ''

  call PetscLogStagePop(ierr)
  
end subroutine Output

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
 
  if (realization%output_option%print_tecplot .or. &
      realization%output_option%print_hdf5) then
    call OutputObservationTecplot(realization)
  endif

end subroutine OutputObservation

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
  use Option_module
  use Field_module
  use Patch_module
  
  use Mphase_module
  use Immis_module
  use THC_module
  use Richards_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXHEADERLENGTH) :: string, string2
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
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.tec'
  else  
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-' // trim(string) // '.tec'
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)
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
      case (IMS_MODE)
        string2 = ImmisGetTecplotHeader(realization,icolumn)
      case (MPH_MODE)
        string2 = MphaseGetTecplotHeader(realization,icolumn)
      case(THC_MODE)
        string2 = THCGetTecplotHeader(realization,icolumn)
      case(RICHARDS_MODE)
       string2 = RichardsGetTecplotHeader(realization,icolumn)
    end select
    string = trim(string) // trim(string2)
    
    ! write transport variables
    if (option%ntrandof > 0) then
      string2 = RTGetTecplotHeader(realization,icolumn)
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
      if (comma_count > 3) then
        write(string2,'(i5)') comma_count+1
        string = trim(string) // ', VARLOCATION=([4-' // &
                 trim(adjustl(string2)) // ']=CELLCENTERED)'
      else
        string = trim(string) // ', VARLOCATION=([4]=CELLCENTERED)'
      endif
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
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)

      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! pressure
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
     ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          ! liquid mole fractions
          do i=1,option%nflowspec
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select
  
      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nflowspec
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select 
#if 0      
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,ZERO_INTEGER)
        call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      endif
#endif    
      
    case default
  
  end select
  
  if (option%ntrandof > 0) then
    if (associated(reaction)) then
      if (reaction%print_pH .and. reaction%h_ion_id > 0) then
        call OutputGetVarFromArray(realization,global_vec,PH,reaction%h_ion_id)
        call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      endif
      do i=1,reaction%ncomp
        if (reaction%primary_species_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,TOTAL_MOLARITY,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
        endif
      enddo
      if (realization%output_option%print_act_coefs) then
        do i=1,reaction%ncomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_ACTIVITY_COEF,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          endif
        enddo
      endif
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
        endif
      enddo
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_RATE,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
        endif
      enddo
      do i=1,reaction%neqsurfcmplxrxn
        if (reaction%surface_site_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX_FREE,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
        endif
      enddo
      do i=1,reaction%neqsurfcmplx
        if (reaction%surface_complex_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX,i)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
        endif
      enddo
    endif
  endif
  
  ! material id
  if (associated(patch%imat)) then
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
  endif

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  if (option%myrank == option%io_rank) close(IUNIT3)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplotBlock(realization)
  endif
  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          X_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization,GAS_PHASE, &
                                              X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE)
          call OutputFluxVelocitiesTecplotBlk(realization,GAS_PHASE, &
                                              Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplotBlk(realization,LIQUID_PHASE, &
                                          Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE, IMS_MODE)
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
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-vel.tec'
  else  
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-vel-' // trim(string) // '.tec'
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity output file: ' // &
                       trim(filename)
    call printMsg(option)
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
             '"vlx [m/' // trim(output_option%tunit) // ']",' // &
             '"vly [m/' // trim(output_option%tunit) // ']",' // &
             '"vlz [m/' // trim(output_option%tunit) // ']"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/' // trim(output_option%tunit) // ']",' // &
               '"vgy [m/' // trim(output_option%tunit) // ']",' // &
               '"vgz [m/' // trim(output_option%tunit) // ']"'
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
      if (associated(patch%imat)) then
        if (option%nphase > 1) then
          string = trim(string) // ', VARLOCATION=([4-10]=CELLCENTERED)'
        else
          string = trim(string) // ', VARLOCATION=([4-7]=CELLCENTERED)'
        endif
      else
        if (option%nphase > 1) then
          string = trim(string) // ', VARLOCATION=([4-9]=CELLCENTERED)'
        else
          string = trim(string) // ', VARLOCATION=([4-6]=CELLCENTERED)'
        endif
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

  if (option%myrank == option%io_rank) close(IUNIT3)
  
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
  
  if (output_option%plot_number < 10) then
    write(string,'("00",i1)') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(string,'("0",i2)') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(string,'(i3)') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(string,'(i4)') output_option%plot_number  
  endif
  
  filename = trim(filename) // '-' // trim(string) // '.tec'
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity flux output file: ' // &
                       trim(filename)
    call printMsg(option)
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
        string = trim(string) // ' vlx [m/' // trim(output_option%tunit) // ']"'
      case(Y_DIRECTION)
        string = trim(string) // ' vly [m/' // trim(output_option%tunit) // ']"'
      case(Z_DIRECTION)
        string = trim(string) // ' vlz [m/' // trim(output_option%tunit) // ']"'
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
        ghosted_id = grid%nL2G(local_id)        
        array(count) = grid%y(ghosted_id)
        if (direction == Y_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dy(ghosted_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
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
        ghosted_id = grid%nL2G(local_id)        
        array(count) = grid%z(ghosted_id)
        if (direction == Z_DIRECTION) &
          array(count) = array(count) + 0.5d0*grid%structured_grid%dz(ghosted_id)
      enddo
    enddo
  enddo
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size,option)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  if (option%myrank == option%io_rank) close(IUNIT3)

  call PetscLogEventEnd(logging%event_output_write_flux_tecplot, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  
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
  use Richards_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  PetscInt :: icolumn
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXHEADERLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch 
  type(reaction_type), pointer :: reaction 
  type(output_option_type), pointer :: output_option
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: value
  Vec :: global_vec
  Vec :: natural_vec
  
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.tec'
  else
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-' // trim(string) // '.tec'    
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // &
                       trim(filename)
    call printMsg(option)                       
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

    icolumn = 3
    
    ! write flow variables
    string2 = ''
    select case(option%iflowmode)
      case (IMS_MODE)
        string2 = ImmisGetTecplotHeader(realization,icolumn)
      case (MPH_MODE)
        string2 = MphaseGetTecplotHeader(realization,icolumn)
      case(THC_MODE)
        string2 = THCGetTecplotHeader(realization,icolumn)
      case(RICHARDS_MODE)
       string2 = RichardsGetTecplotHeader(realization,icolumn)
    end select
    string = trim(string) // trim(string2)
    
    ! write transport variables
    if (option%ntrandof > 0) then
      string2 = RTGetTecplotHeader(realization,icolumn)
      string = trim(string) // trim(string2)
    endif

    ! write material ids
    if (associated(patch%imat)) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Material_ID"'')') icolumn
      string = trim(string) // trim(string2)
    endif

    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i5,'', J='',i5, &
                   &'', K='',i5)') &
                   option%time/output_option%tconv,grid%structured_grid%nx, &
                   grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ', DATAPACKING=POINT'
    write(IUNIT3,'(a)') trim(string)
  endif
  
1000 format(es13.6,1x)
1001 format(i11,1x)
1009 format('')

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    write(IUNIT3,1000,advance='no') grid%x(ghosted_id)
    write(IUNIT3,1000,advance='no') grid%y(ghosted_id)
    write(IUNIT3,1000,advance='no') grid%z(ghosted_id)

    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)

        ! temperature
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,TEMPERATURE, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select

        ! pressure
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,PRESSURE, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select

        ! phase
        select case(option%iflowmode)
          case(MPH_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,PHASE, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1001,advance='no') int(value)
        end select

        ! liquid saturation
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,LIQUID_SATURATION, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select

        ! gas saturation
        select case(option%iflowmode)
          case(MPH_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,GAS_SATURATION, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select
      
        ! liquid energy
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,LIQUID_ENERGY, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select
      
       ! gas energy
        select case(option%iflowmode)
          case(MPH_MODE,IMS_MODE)
            value = RealizGetDatasetValueAtCell(realization,GAS_ENERGY, &
                                                ZERO_INTEGER,ghosted_id)
            write(IUNIT3,1000,advance='no') value
        end select

        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE)
            ! liquid mole fractions
            do i=1,option%nflowspec
              value = RealizGetDatasetValueAtCell(realization,LIQUID_MOLE_FRACTION, &
                                                  i,ghosted_id)
              write(IUNIT3,1000,advance='no') value
            enddo
        end select
    
        select case(option%iflowmode)
          case(MPH_MODE)
            ! gas mole fractions
            do i=1,option%nflowspec
              value = RealizGetDatasetValueAtCell(realization,GAS_MOLE_FRACTION, &
                                                  i,ghosted_id)
              write(IUNIT3,1000,advance='no') value
            enddo
        end select 

      case default
    
    end select
    
    if (option%ntrandof > 0) then
      if (associated(reaction)) then
        if (reaction%print_pH .and. reaction%h_ion_id > 0) then
          value = RealizGetDatasetValueAtCell(realization,PH, &
                                              reaction%h_ion_id,ghosted_id)
          write(IUNIT3,1000,advance='no') value
        endif
        do i=1,reaction%ncomp
          if (reaction%primary_species_print(i)) then
            value = RealizGetDatasetValueAtCell(realization,TOTAL_MOLARITY, &
                                                i,ghosted_id)
            write(IUNIT3,1000,advance='no') value
          endif
        enddo
        if (realization%output_option%print_act_coefs) then
          do i=1,reaction%ncomp
            if (reaction%primary_species_print(i)) then
              value = RealizGetDatasetValueAtCell(realization,PRIMARY_ACTIVITY_COEF, &
                                                  i,ghosted_id)
              write(IUNIT3,1000,advance='no') value
            endif
          enddo
        endif
        do i=1,reaction%nkinmnrl
          if (reaction%kinmnrl_print(i)) then
            value = RealizGetDatasetValueAtCell(realization,MINERAL_VOLUME_FRACTION, &
                                                i,ghosted_id)
            write(IUNIT3,1000,advance='no') value
          endif
        enddo
        do i=1,reaction%nkinmnrl
          if (reaction%kinmnrl_print(i)) then
            value = RealizGetDatasetValueAtCell(realization,MINERAL_RATE, &
                                                i,ghosted_id)
            write(IUNIT3,1000,advance='no') value
          endif
        enddo
        do i=1,reaction%neqsurfcmplxrxn
          if (reaction%surface_site_print(i)) then
            value = RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX_FREE, &
                                                i,ghosted_id)
            write(IUNIT3,1000,advance='no') value
          endif
        enddo
        do i=1,reaction%neqsurfcmplx
          if (reaction%surface_complex_print(i)) then
            value = RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX, &
                                                i,ghosted_id)
            write(IUNIT3,1000,advance='no') value
          endif
        enddo
      endif
    endif
    
    ! material id
    if (associated(patch%imat)) then
      value = RealizGetDatasetValueAtCell(realization,MATERIAL_ID, &
                                              ZERO_INTEGER,ghosted_id)
      write(IUNIT3,1001,advance='no') int(value)
    endif
    write(IUNIT3,1009) 

  enddo
  
  
  if (option%myrank == option%io_rank) close(IUNIT3)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplotPoint(realization)
  endif
  
end subroutine OutputTecplotPoint

! ************************************************************************** !
!
! OutputVelocitiesTecplotPoint: Print velocities to Tecplot file in BLOCK format
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
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-vel.tec'
  else  
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-vel-' // trim(string) // '.tec'
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot velocity output file: ' // &
                       trim(filename)
    call printMsg(option)                       
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
             '"vlx [m/' // trim(output_option%tunit) // ']",' // &
             '"vly [m/' // trim(output_option%tunit) // ']",' // &
             '"vlz [m/' // trim(output_option%tunit) // ']"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/' // trim(output_option%tunit) // ']",' // &
               '"vgy [m/' // trim(output_option%tunit) // ']",' // &
               '"vgz [m/' // trim(output_option%tunit) // ']"'
    endif
        if (associated(patch%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i5,'', J='',i5, &
                 &'', K='',i5,'','')') &
                 option%time/output_option%tconv, &
                 grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=POINT'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! currently supported for only liquid phase
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vx,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vy,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec_vz,GLOBAL, &
                                  option)  
  
  call GetCellCenteredVelocities(realization,global_vec_vx,LIQUID_PHASE,X_DIRECTION)
  call GetCellCenteredVelocities(realization,global_vec_vy,LIQUID_PHASE,Y_DIRECTION)
  call GetCellCenteredVelocities(realization,global_vec_vz,LIQUID_PHASE,Z_DIRECTION)

  call GridVecGetArrayF90(grid,global_vec_vx,vec_ptr_vx,ierr)
  call GridVecGetArrayF90(grid,global_vec_vy,vec_ptr_vy,ierr)
  call GridVecGetArrayF90(grid,global_vec_vz,vec_ptr_vz,ierr)

  ! write points
1000 format(es13.6,1x)
1001 format(i11,1x)
1002 format(3(es13.6,1x))
1009 format('')

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)  ! local and ghosted are same for non-parallel
    write(IUNIT3,1000,advance='no') grid%x(ghosted_id)
    write(IUNIT3,1000,advance='no') grid%y(ghosted_id)
    write(IUNIT3,1000,advance='no') grid%z(ghosted_id)
    
    write(IUNIT3,1000,advance='no') vec_ptr_vx(ghosted_id)
    write(IUNIT3,1000,advance='no') vec_ptr_vy(ghosted_id)
    write(IUNIT3,1000,advance='no') vec_ptr_vz(ghosted_id)

    ! material id
    if (associated(patch%imat)) then
      value = RealizGetDatasetValueAtCell(realization,MATERIAL_ID, &
                                          ZERO_INTEGER,ghosted_id)
      write(IUNIT3,1001,advance='no') int(value)
    endif
  
    write(IUNIT3,1009)
    
  enddo
  
  call GridVecRestoreArrayF90(grid,global_vec_vx,vec_ptr_vx,ierr)
  call GridVecRestoreArrayF90(grid,global_vec_vy,vec_ptr_vy,ierr)
  call GridVecRestoreArrayF90(grid,global_vec_vz,vec_ptr_vz,ierr)
  
  call VecDestroy(global_vec_vx,ierr)
  call VecDestroy(global_vec_vy,ierr)
  call VecDestroy(global_vec_vz,ierr)

  if (option%myrank == option%io_rank) close(IUNIT3)
  
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

  call PetscLogEventBegin(logging%event_output_vec_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 

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
    if (associated(patch%imat)) then
      string = trim(string) // ',"Material_ID"'
    endif
    write(fid,'(a)') trim(string)
  
    ! write zone header
    if (realization%discretization%itype == STRUCTURED_GRID) then
      write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i5,'', J='',i5, &
                   &'', K='',i5,'','')') &
                   option%time/realization%output_option%tconv, &
                   grid%structured_grid%nx+1,grid%structured_grid%ny+1,grid%structured_grid%nz+1
      string = trim(string) // ' DATAPACKING=BLOCK'
                                                ! 4=dataset name, 5=material_id
      if (associated(patch%imat)) then
        string = trim(string) // ', VARLOCATION=([4-5]=CELLCENTERED)'
      else
        string = trim(string) // ', VARLOCATION=([4]=CELLCENTERED)'
      endif
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

  if (realization%discretization%itype == STRUCTURED_GRID) then
    call WriteTecplotStructuredGrid(fid,realization)
  else  
    call GetCoordinates(grid,global_vec,X_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Y_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

    call GetCoordinates(grid,global_vec,Z_COORDINATE)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)
  endif    

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

1000 format(es13.6,1x)
1001 format(10(es13.6,1x))
  
  call PetscLogEventBegin(logging%event_output_str_grid_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                              
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
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc, recv_size
  PetscInt :: max_local_size, local_size
  PetscInt :: istart, iend, num_in_array
  PetscInt :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)

1000 format(es13.6)
1001 format(10(es13.6,1x))
!1000 format(es16.9)
!1001 format(10(es16.9,1x))
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call PetscLogEventBegin(logging%event_output_write_tecplot, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

  ! maximum number of initial messages  
#define HANDSHAKE  
  max_proc = option%io_handshake_buffer_size
  max_proc_prefetch = option%io_handshake_buffer_size / 10

  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
    local_size = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                         option%mycomm,ierr)
      max_local_size_saved = max_local_size
      write(option%io_buffer,'("max_local_size_saved: ",i9)') max_local_size
      call printMsg(option)
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
  if (option%myrank == option%io_rank) then
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
        write(fid,1001) real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size-iend) = real_data(iend+1:local_size)
      num_in_array = local_size-iend
    endif
    do iproc=1,option%mycommsize-1
#ifdef HANDSHAKE    
      if (option%io_handshake_buffer_size > 0 .and. &
          iproc+max_proc_prefetch >= max_proc) then
        max_proc = max_proc + option%io_handshake_buffer_size
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      endif
#endif      
      call MPI_Probe(iproc,MPI_ANY_TAG,option%mycomm,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data_recv,recv_size,MPI_INTEGER,iproc, &
                      MPI_ANY_TAG,option%mycomm,status,ierr)
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
                      MPI_ANY_TAG,option%mycomm,status,ierr)
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
      call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                     ierr)
    endif
#endif      
    ! Print the remaining values, if they exist
    if (datatype == 0) then
      if (num_in_array > 0) &
        write(fid,'(10(i3,x))') integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(fid,1001) real_data(1:num_in_array)
    endif
  else
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        if (option%myrank < max_proc) exit
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      enddo
    endif
#endif    
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size,MPI_INTEGER,option%io_rank,local_size, &
                    option%mycomm,ierr)
    else
      call MPI_Send(real_data,local_size,MPI_DOUBLE_PRECISION,option%io_rank,local_size, &
                    option%mycomm,ierr)
    endif
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
        if (max_proc < 0) exit
      enddo
    endif
#endif    
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
  PetscTruth, save :: open_file = PETSC_FALSE
  PetscInt :: local_id

  call PetscLogEventBegin(logging%event_output_observation, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  if (observation_first) then
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
    enddo
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
    if (observation_first) then
      open(unit=fid,file=filename,action="write",status="replace")
      ! write header
      ! write title
      write(fid,'(a)',advance="no") '"Time[' // trim(output_option%tunit) // ']"'
      observation => patch%observation%first
      do 
        if (.not.associated(observation)) exit
        select case(observation%itype)
          case(OBSERVATION_SCALAR)
            if (associated(observation%region%coordinates)) then
              call WriteObservationHeaderForCoord(fid,realization, &
                                                   observation%region, &
                                                   observation%print_velocities)
            else
              do icell=1,observation%region%num_cells
                call WriteObservationHeaderForCell(fid,realization, &
                                                    observation%region,icell, &
                                                    observation%print_velocities)
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
            if (associated(observation%region%coordinates)) then
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
  
  call PetscLogEventEnd(logging%event_output_observation, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
      
end subroutine OutputObservationTecplot

! ************************************************************************** !
!
! WriteObservationHeaderForCell: Print a header for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteObservationHeaderForCell(fid,realization,region,icell, &
                                          print_velocities)

  use Realization_module
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module
  use Region_module
  use Reaction_Aux_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscInt :: icell
  PetscTruth :: print_velocities
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: cell_string
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction  
  
  patch => realization%patch
  option => realization%option
  field => realization%field
  grid => patch%grid
  
  write(cell_string,*) grid%nL2A(region%cell_ids(icell))
  cell_string = trim(region%name) // ' ' //adjustl(cell_string)

  select case(option%iflowmode)
    case (IMS_MODE)
!      string = ',"X [m] '// trim(cell_string) // '",' // &
!               '"Y [m] '// trim(cell_string) // '",' // &
!               '"Z [m] '// trim(cell_string) // '",' // &
      string = '"T [C] '// trim(cell_string) // '",' // &
               '"P [Pa] '// trim(cell_string) // '",' // &
               '"sl '// trim(cell_string) // '",' // &
               '"sg '// trim(cell_string) // '",' // &
               '"Ul '// trim(cell_string) // '",' // &
               '"Ug '// trim(cell_string) // '",'
    case (MPH_MODE)
!      string = ',"X [m] '// trim(cell_string) // '",' // &
!               '"Y [m] '// trim(cell_string) // '",' // &
!               '"Z [m] '// trim(cell_string) // '",' // &
      string = '"T [C] '// trim(cell_string) // '",' // &
               '"P [Pa] '// trim(cell_string) // '",' // &
               '"sl '// trim(cell_string) // '",' // &
               '"sg '// trim(cell_string) // '",' // &
               '"Ul '// trim(cell_string) // '",' // &
               '"Ug '// trim(cell_string) // '",'
      do i=1,option%nflowspec
        write(string2,'(''"Xl('',i2,'') '// trim(cell_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,option%nflowspec
        write(string2,'(''"Xg('',i2,'') '// trim(cell_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      string = trim(string) // ',"Phase '// trim(cell_string) // '"'
    case(THC_MODE,RICHARDS_MODE)
      if (option%iflowmode == THC_MODE) then
!        string = ',"X [m] '// trim(cell_string) // '",' // &
!                 '"Y [m] '// trim(cell_string) // '",' // &
!                 '"Z [m] '// trim(cell_string) // '",' // &
        string = '"T [C] '// trim(cell_string) // '",' // &
                 '"P [Pa] '// trim(cell_string) // '",' // &
                 '"sl '// trim(cell_string) // '",' // &
                 '"Ul '// trim(cell_string) // '"' 
      else
!        string = ',"X [m] '// trim(cell_string) // '",' // &
!                 '"Y [m] '// trim(cell_string) // '",' // &
!                 '"Z [m] '// trim(cell_string) // '",' // &
        string = ',"P [Pa] '// trim(cell_string) // '",' // &
                 '"sl '// trim(cell_string) // '"'
      endif
      if (option%iflowmode == THC_MODE) then
        do i=1,option%nflowspec
          write(string2,'('',"Xl('',i2,'') '// trim(cell_string) // '"'')') i
          string = trim(string) // trim(string2)
        enddo
      endif
    case default
      string = ''
  end select
  write(fid,'(a)',advance="no") trim(string)

  ! reactive transport
  if (option%ntrandof > 0) then
 
    reaction => realization%reaction
    if ((reaction%print_pH) .and. &
        reaction%h_ion_id > 0) then
      write(fid,'('',"pH '',a,''"'')',advance="no") trim(cell_string)
    endif
    
    do i=1,option%ntrandof
      if (reaction%primary_species_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%primary_species_names(i)), trim(cell_string)
      endif
    enddo
    
    if (realization%output_option%print_act_coefs) then
      do i=1,option%ntrandof
        if (reaction%primary_species_print(i)) then
          write(fid,'('',"'',a,''_gam '',a,''"'')',advance="no") &
            trim(reaction%primary_species_names(i)), trim(cell_string)
        endif
      enddo
    endif
    
    do i=1,reaction%nkinmnrl
      if (reaction%kinmnrl_print(i)) then
        write(fid,'('',"'',a,''_vf '',a,''"'')',advance="no") &
          trim(reaction%kinmnrl_names(i)), trim(cell_string)
      endif
    enddo
    
    do i=1,reaction%nkinmnrl
      if (reaction%kinmnrl_print(i)) then
        write(fid,'('',"'',a,''_rt '',a,''"'')',advance="no") &
          trim(reaction%kinmnrl_names(i)), trim(cell_string)  
      endif
    enddo
    
    do i=1,realization%reaction%neqsurfcmplxrxn
      if (reaction%surface_site_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%surface_site_names(i)), trim(cell_string)
      endif
    enddo
    
    do i=1,realization%reaction%neqsurfcmplx
      if (reaction%surface_complex_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%surface_complex_names(i)), trim(cell_string)
      endif
    enddo

  endif

  if (print_velocities) then 
    string = ',"vlx [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"' // &
             ',"vly [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"' // &
             ',"vlz [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"'
    write(fid,'(a)',advance="no") trim(string)
  endif

end subroutine WriteObservationHeaderForCell

! ************************************************************************** !
!
! WriteObservationHeaderForCoord: Print a header for data at a coordinate
! author: Glenn Hammond
! date: 04/11/08
!
! ************************************************************************** !  
subroutine WriteObservationHeaderForCoord(fid,realization,region, &
                                           print_velocities)

  use Realization_module
  use Option_module
  use Patch_module
  use Region_module
  use Reaction_Aux_module
  
  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  type(region_type) :: region
  PetscTruth :: print_velocities
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXSTRINGLENGTH) :: coordinate_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction  
  
  patch => realization%patch
  option => realization%option
  
!  write(cell_string,*) grid%nL2A(region%cell_ids(icell))
!  cell_string = trim(region%name) // ' ' //adjustl(cell_string)
  cell_string = trim(region%name)
  
  110 format(1pg12.4)
  write(x_string,110) region%coordinates(ONE_INTEGER)%x
  write(y_string,110) region%coordinates(ONE_INTEGER)%y
  write(z_string,110) region%coordinates(ONE_INTEGER)%z
  cell_string = trim(cell_string) // ' ' // trim(adjustl(x_string)) // ' ' // &
                   trim(adjustl(y_string)) // ' ' // trim(adjustl(z_string))

  select case(option%iflowmode)
    case (IMS_MODE)
!      string = ',"X [m] '// trim(cell_string) // '",' // &
!               '"Y [m] '// trim(cell_string) // '",' // &
!               '"Z [m] '// trim(cell_string) // '",' // &
       string = ',"T [C] '// trim(cell_string) // '",' // &
               '"P [Pa] '// trim(cell_string) // '",' // &
               '"sl '// trim(cell_string) // '",' // &
               '"sg '// trim(cell_string) // '",' // &
               '"Ul '// trim(cell_string) // '",' // &
               '"Ug '// trim(cell_string) // '",'
    case (MPH_MODE)
!      string = ',"X [m] '// trim(cell_string) // '",' // &
!               '"Y [m] '// trim(cell_string) // '",' // &
!               '"Z [m] '// trim(cell_string) // '",' // &
       string = ',"T [C] '// trim(cell_string) // '",' // &
               '"P [Pa] '// trim(cell_string) // '",' // &
               '"sl '// trim(cell_string) // '",' // &
               '"sg '// trim(cell_string) // '",' // &
               '"Ul '// trim(cell_string) // '",' // &
               '"Ug '// trim(cell_string) // '",'
      do i=1,option%nflowspec
        write(string2,'(''"Xl('',i2,'') '// trim(cell_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,option%nflowspec
        write(string2,'(''"Xg('',i2,'') '// trim(cell_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      string = trim(string) // ',"Phase '// trim(cell_string) // '"'
    case(THC_MODE,RICHARDS_MODE)
      if (option%iflowmode == THC_MODE) then
!        string = ',"X [m] '// trim(cell_string) // '",' // &
!                 '"Y [m] '// trim(cell_string) // '",' // &
!                 '"Z [m] '// trim(cell_string) // '",' // &
        string = ',"T [C] '// trim(cell_string) // '",' // &
                 '"P [Pa] '// trim(cell_string) // '",' // &
                 '"sl '// trim(cell_string) // '",' // &
                 '"Ul '// trim(cell_string) // '"' 
      else
!        string = ',"X [m] '// trim(cell_string) // '",' // &
!                 '"Y [m] '// trim(cell_string) // '",' // &
!                 '"Z [m] '// trim(cell_string) // '",' // &
        string = ',"P [Pa] '// trim(cell_string) // '",' // &
                 '"sl '// trim(cell_string) // '"'
      endif
      if (option%iflowmode == THC_MODE) then
        do i=1,option%nflowspec
          write(string2,'('',"Xl('',i2,'') '// trim(cell_string) // '"'')') i
          string = trim(string) // trim(string2)
        enddo
      endif
    case default
      string = ''
  end select
  write(fid,'(a)',advance="no") trim(string)

  ! reactive transport
  if (option%ntrandof > 0) then

    reaction => realization%reaction
    if ((reaction%print_pH) .and. &
        reaction%h_ion_id > 0) then
      write(fid,'('',"pH '',a,''"'')',advance="no") trim(cell_string)
    endif
    
    do i=1,option%ntrandof
      if (reaction%primary_species_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%primary_species_names(i)), trim(cell_string)
      endif
    enddo
    
    if (realization%output_option%print_act_coefs) then
      do i=1,option%ntrandof
        if (reaction%primary_species_print(i)) then
          write(fid,'('',"'',a,''_gam '',a,''"'')',advance="no") &
            trim(reaction%primary_species_names(i)), trim(cell_string)
        endif
      enddo
    endif
    
    do i=1,reaction%nkinmnrl
      if (reaction%kinmnrl_print(i)) then
        write(fid,'('',"'',a,''_vf '',a,''"'')',advance="no") &
          trim(reaction%kinmnrl_names(i)), trim(cell_string)
      endif
    enddo
    
    do i=1,reaction%nkinmnrl
      if (reaction%kinmnrl_print(i)) then
        write(fid,'('',"'',a,''_rt '',a,''"'')',advance="no") &
          trim(reaction%kinmnrl_names(i)), trim(cell_string)  
      endif
    enddo
    
    do i=1,realization%reaction%neqsurfcmplxrxn
      if (reaction%surface_site_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%surface_site_names(i)), trim(cell_string)
      endif
    enddo
    
    do i=1,realization%reaction%neqsurfcmplx
      if (reaction%surface_complex_print(i)) then
        write(fid,'('',"'',a,'' '',a,''"'')',advance="no") &
          trim(reaction%surface_complex_names(i)), trim(cell_string)
      endif
    enddo

  endif

  if (print_velocities) then 
    string = ',"vlx [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"' // &
             ',"vly [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"' // &
             ',"vlz [m/'//trim(realization%output_option%tunit)//'] '// &
             trim(cell_string) // '"'
    write(fid,'(a)',advance="no") trim(string)
  endif

end subroutine WriteObservationHeaderForCoord

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
    case (MPH_MODE)
    case (IMS_MODE)
    case(THC_MODE)
    case(RICHARDS_MODE)
      string = ',"Darcy flux ' // trim(coupler_name) // &
               ' [m^3/' // trim(realization%output_option%tunit) // ']"'
    case default
  end select
  write(fid,'(a)',advance="no") trim(string)

  if (associated(reaction)) then
    do i=1, reaction%ncomp 
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
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

100 format(es14.6)
101 format("i")
110 format(es14.6)
111 format("i")

  ghosted_id = grid%nL2G(local_id)
  ! write out coorindates
  !write(fid,110,advance="no") grid%x(ghosted_id)
  !write(fid,110,advance="no") grid%y(ghosted_id)
  !write(fid,110,advance="no") grid%z(ghosted_id)
  
  ! temperature
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,TEMPERATURE,ZERO_INTEGER,ghosted_id)
  end select

  ! pressure
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,PRESSURE,ZERO_INTEGER,ghosted_id)
  end select

  ! liquid saturation
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,LIQUID_SATURATION,ZERO_INTEGER,ghosted_id)
  end select

 ! gas saturation
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,GAS_SATURATION,ZERO_INTEGER,ghosted_id)
  end select

  ! liquid energy
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,LIQUID_ENERGY,ZERO_INTEGER,ghosted_id)
  end select

  ! gas energy
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        RealizGetDatasetValueAtCell(realization,GAS_ENERGY,ZERO_INTEGER,ghosted_id)
  end select

  ! liquid mole fractions
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      do i=1,option%nflowspec
        write(fid,110,advance="no") &
          RealizGetDatasetValueAtCell(realization,LIQUID_MOLE_FRACTION,i,ghosted_id)
      enddo
  end select

  ! gas mole fractions
  select case(option%iflowmode)
    case(MPH_MODE)
      do i=1,option%nflowspec
        write(fid,110,advance="no") &
          RealizGetDatasetValueAtCell(realization,GAS_MOLE_FRACTION,i,ghosted_id)
      enddo
  end select 

  ! phase
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,111,advance="no") &
        int(RealizGetDatasetValueAtCell(realization,PHASE,ZERO_INTEGER,ghosted_id))
  end select

  if (option%ntrandof > 0) then
    reaction => realization%reaction
    ghosted_id = grid%nL2G(local_id)
    if (associated(reaction)) then
      if (reaction%print_pH .and. reaction%h_ion_id > 0) then
        write(fid,110,advance="no") &
          RealizGetDatasetValueAtCell(realization,PH,reaction%h_ion_id,ghosted_id)
      endif
      do i=1,reaction%ncomp
        if (reaction%primary_species_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,TOTAL_MOLARITY,i,ghosted_id)
        endif
      enddo
      if (realization%output_option%print_act_coefs) then
        do i=1,reaction%ncomp
          if (reaction%primary_species_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,PRIMARY_ACTIVITY_COEF,i,ghosted_id)
          endif
        enddo
      endif
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,MINERAL_VOLUME_FRACTION,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,MINERAL_RATE,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%neqsurfcmplxrxn
        if (reaction%surface_site_print(i)) then
           write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX_FREE,i,ghosted_id)
        endif
      enddo
      do i=1,reaction%neqsurfcmplx
        if (reaction%surface_complex_print(i)) then
          write(fid,110,advance="no") &
            RealizGetDatasetValueAtCell(realization,SURFACE_CMPLX,i,ghosted_id)
        endif
      enddo
    endif
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
  PetscInt :: ghosted_ids(8)
  PetscInt :: count
  PetscInt :: i, j, k
  PetscInt :: istart, iend, jstart, jend, kstart, kend
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

100 format(es14.6)
!100 format(es16.9)
101 format("i")
110 format(es14.6)
!110 format(',',es16.9)
111 format("i")


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
  
  ! temperature
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,TEMPERATURE,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! pressure
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,PRESSURE,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! liquid saturation
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,LIQUID_SATURATION,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! gas saturation
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE)
      ! gas saturation
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,GAS_SATURATION,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! liquid energy
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,LIQUID_ENERGY,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! gas energy
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE)
      write(fid,110,advance="no") &
        OutputGetVarFromArrayAtCoord(realization,GAS_ENERGY,ZERO_INTEGER, &
                                     region%coordinates(ONE_INTEGER)%x, &
                                     region%coordinates(ONE_INTEGER)%y, &
                                     region%coordinates(ONE_INTEGER)%z, &
                                     count,ghosted_ids)
  end select

  ! liquid mole fraction
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
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
    case(MPH_MODE)
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
    case(MPH_MODE,THC_MODE,IMS_MODE)
      write(fid,111,advance="no") &
        int(OutputGetVarFromArrayAtCoord(realization,PHASE,ZERO_INTEGER, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids))
  end select

  if (option%ntrandof > 0) then
    reaction => realization%reaction
    if (associated(reaction)) then
      if (reaction%print_pH .and. reaction%h_ion_id > 0) then
        write(fid,110,advance="no") &
          OutputGetVarFromArrayAtCoord(realization,PH,reaction%h_ion_id, &
                                       region%coordinates(ONE_INTEGER)%x, &
                                       region%coordinates(ONE_INTEGER)%y, &
                                       region%coordinates(ONE_INTEGER)%z, &
                                       count,ghosted_ids)
      endif
      do i=1,reaction%ncomp
        if (reaction%primary_species_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,TOTAL_MOLARITY,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      if (realization%output_option%print_act_coefs) then
        do i=1,reaction%ncomp
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
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,MINERAL_VOLUME_FRACTION,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,MINERAL_RATE,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%neqsurfcmplxrxn
        if (reaction%surface_site_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,SURFACE_CMPLX_FREE,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
      do i=1,reaction%neqsurfcmplx
        if (reaction%surface_complex_print(i)) then
          write(fid,110,advance="no") &
            OutputGetVarFromArrayAtCoord(realization,SURFACE_CMPLX,i, &
                                         region%coordinates(ONE_INTEGER)%x, &
                                         region%coordinates(ONE_INTEGER)%y, &
                                         region%coordinates(ONE_INTEGER)%z, &
                                         count,ghosted_ids)
        endif
      enddo
    endif
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
  PetscReal :: sum_volumetric_flux(realization%option%nphase)
  PetscReal :: sum_volumetric_flux_global(realization%option%nphase)
  PetscReal :: sum_solute_flux(realization%reaction%ncomp)
  PetscReal :: sum_solute_flux_global(realization%reaction%ncomp)
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  
  option => realization%option
  reaction => realization%reaction

100 format(es14.6)
!100 format(es16.9)
101 format("i")
110 format(es14.6)
!110 format(',',es16.9)
111 format("i")
 
  iphase = 1

  ! sum up fluxes across region
  if (associated(connection_set)) then
    offset = connection_set%offset
    select case(option%iflowmode)
      case(MPH_MODE,THC_MODE,IMS_MODE)
      case(RICHARDS_MODE)
        sum_volumetric_flux = 0.d0
        if (associated(connection_set)) then
          do iconn = 1, connection_set%num_connections
            sum_volumetric_flux(:) = sum_volumetric_flux(:) + &
                                  patch%boundary_velocities(iphase,offset+iconn)* &
                                  connection_set%area(iconn)
          enddo
        endif
        call MPI_Reduce(sum_volumetric_flux,sum_volumetric_flux_global, &
                        option%nphase,MPI_DOUBLE_PRECISION,MPI_SUM, &
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
      call MPI_Reduce(sum_solute_flux,sum_solute_flux_global, &
                      reaction%ncomp,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
      if (option%myrank == option%io_rank) then
        do i = 1, reaction%ncomp
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
                 dabs(cur_connection_set%dist(direction,iconn))
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + area
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
                 dabs(cur_connection_set%dist(direction,iconn))
          sum_velocity(direction) = sum_velocity(direction) + &
                                    patch%internal_velocities(iphase,sum_connection)* &
                                    area
          sum_area(direction) = sum_area(direction) + area
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
  
  use Mphase_module
  use Immis_module
  use THC_module
  use Richards_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string, string2
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
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.vtk'
  else
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-' // trim(string) // '.vtk'    
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write vtk output file: ' // trim(filename)
    call printMsg(option)    
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    write(IUNIT3,'(''# vtk DataFile Version 2.0'')')
    ! write title
    write(IUNIT3,'(''PFLOTRAN output'')')
    write(IUNIT3,'(''ASCII'')')
    write(IUNIT3,'(''DATASET UNSTRUCTURED_GRID'')')
  endif

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)  
  call DiscretizationCreateVector(discretization,ONEDOF,natural_vec,NATURAL, &
                                  option)  

  ! write out coordinates
  call WriteVTKGrid(IUNIT3,realization)

  write(IUNIT3,'(''CELL_DATA'',i8)') grid%nmax

  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)

      ! temperature
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          word = 'Temperature'
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select

      ! pressure
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          word = 'Pressure'
          call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select

      ! phase
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          word = 'Phase'
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_INTEGER)
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          word = 'Liquid_Saturation'
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select

      ! gas saturation
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          word = 'Gas_Saturation'
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select
    
      ! liquid energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          word = 'Liquid_Energy'
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select
    
     ! gas energy
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          word = 'Gas_Energy'
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
      end select

      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE)
          ! liquid mole fractions
          do i=1,option%nflowspec
            write(word,'(''Xl('',i2,'')'')') i
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
          enddo
      end select
  
      select case(option%iflowmode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nflowspec
            write(word,'(''Xg('',i2,'')'')') i
            call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
            call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
            call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_REAL)
          enddo
      end select 
      
    case default
  
  end select
  
  if (option%ntrandof > 0) then
    if (associated(reaction)) then
      do i=1,reaction%ncomp
        call OutputGetVarFromArray(realization,global_vec,TOTAL_MOLARITY,i)
        call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
        call WriteVTKDataSetFromVec(IUNIT3,realization,reaction%primary_species_names(i), &
                                    natural_vec,VTK_REAL)
      enddo
      do i=1,reaction%nkinmnrl
        call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,i)
        call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
        call WriteVTKDataSetFromVec(IUNIT3,realization,reaction%kinmnrl_names(i), &
                                    natural_vec,VTK_REAL)
      enddo
    endif
  endif
  
  ! material id
  if (associated(patch%imat)) then
    word = 'Material_ID'
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteVTKDataSetFromVec(IUNIT3,realization,word,natural_vec,VTK_INTEGER)
  endif

  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  if (option%myrank == option%io_rank) close(IUNIT3)

#if 0  
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
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesVTK(realization,LIQUID_PHASE, &
                                          Y_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesVTK(realization,LIQUID_PHASE, &
                                          Z_DIRECTION)
      select case(option%iflowmode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesVTK(realization,GAS_PHASE, &
                                              Z_DIRECTION)
      end select
    endif
  endif
#endif
      
end subroutine OutputVTK

#if 0
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
    filename = trim(output_option%plot_name) // '-vel.tec'
  else  
    if (output_option%plot_number < 10) then
      write(string,'("00",i1)') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2)') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3)') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4)') output_option%plot_number  
    endif
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-vel-' // trim(string) // '.tec'
  endif
  
  if (option%myrank == option%io_rank) then
   option%io_buffer = '--> write tecplot velocity output file: ' // &
                      trim(filename)
    call printMsg(option)                      
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
             '"vlx [m/' // trim(output_option%tunit) // ']",' // &
             '"vly [m/' // trim(output_option%tunit) // ']",' // &
             '"vlz [m/' // trim(output_option%tunit) // ']"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/' // trim(output_option%tunit) // ']",' // &
               '"vgy [m/' // trim(output_option%tunit) // ']",' // &
               '"vgz [m/' // trim(output_option%tunit) // ']"'
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
      if (associated(patch%imat)) then
        if (option%nphase > 1) then
          string = trim(string) // ', VARLOCATION=([4-10]=CELLCENTERED)'
        else
          string = trim(string) // ', VARLOCATION=([4-7]=CELLCENTERED)'
        endif
      else
        if (option%nphase > 1) then
          string = trim(string) // ', VARLOCATION=([4-9]=CELLCENTERED)'
        else
          string = trim(string) // ', VARLOCATION=([4-6]=CELLCENTERED)'
        endif
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

  ! write out coordinates
  call WriteVTKGrid(IUNIT3,realization)

  word = 'velx
  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)

  if (option%nphase > 1) then
    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_REAL)
  endif

  ! material id
  if (associated(patch%imat)) then
    call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,VTK_INTEGER)
  endif
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if (option%myrank == option%io_rank) close(IUNIT3)
  
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
  
  call PetscLogEventBegin(logging%event_output_grid_vtk, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                              
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  if (realization%discretization%itype == STRUCTURED_GRID) then

    nx = grid%structured_grid%nx
    ny = grid%structured_grid%ny
    nz = grid%structured_grid%nz
  
    nxp1 = nx+1
    nyp1 = ny+1
    nyp1 = nz+1
  
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

  call PetscLogEventEnd(logging%event_output_grid_vtk, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                            
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
  PetscMPIInt :: iproc, recv_size
  PetscInt :: max_local_size, local_size
  PetscInt :: istart, iend, num_in_array
  PetscInt :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)

1001 format(10(es13.6,1x))
1002 format(i3)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  call PetscLogEventBegin(logging%event_output_write_vtk, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

  ! maximum number of initial messages  
#define HANDSHAKE  
  max_proc = option%io_handshake_buffer_size
  max_proc_prefetch = option%io_handshake_buffer_size / 10

  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
    local_size = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                         option%mycomm,ierr)
      max_local_size_saved = max_local_size
      if (OptionPrintToScreen(option)) print *, 'max_local_size_saved: ', max_local_size
    endif
    max_local_size = max_local_size_saved
    local_size = grid%nlmax
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == VTK_INTEGER) then
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
        if (iend+10 > local_size) exit
        iend = istart+9
        write(fid,1002) integer_data(istart:iend)
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
        write(fid,1001) real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size-iend) = real_data(iend+1:local_size)
      num_in_array = local_size-iend
    endif
    do iproc=1,option%mycommsize-1
#ifdef HANDSHAKE    
      if (option%io_handshake_buffer_size > 0 .and. &
          iproc+max_proc_prefetch >= max_proc) then
        max_proc = max_proc + option%io_handshake_buffer_size
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      endif
#endif      
      call MPI_Probe(iproc,MPI_ANY_TAG,option%mycomm,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data_recv,recv_size,MPI_INTEGER,iproc, &
                      MPI_ANY_TAG,option%mycomm,status,ierr)
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
          write(fid,1002) integer_data(istart:iend)
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size,MPI_DOUBLE_PRECISION,iproc, &
                      MPI_ANY_TAG,option%mycomm,status,ierr)
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
      call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                     ierr)
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
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      enddo
    endif
#endif    
    if (datatype == VTK_INTEGER) then
      call MPI_Send(integer_data,local_size,MPI_INTEGER,option%io_rank, &
                    local_size,option%mycomm,ierr)
    else
      call MPI_Send(real_data,local_size,MPI_DOUBLE_PRECISION,option%io_rank, &
                    local_size,option%mycomm,ierr)
    endif
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        call MPI_Bcast(max_proc,1,MPI_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
        if (max_proc < 0) exit
      enddo
    endif
#endif    
  endif
      
  if (datatype == VTK_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

  call PetscLogEventEnd(logging%event_output_write_vtk, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

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
  
  use AMR_Grid_Module
 
#ifndef USE_HDF5
  implicit none
  
  type(realization_type) :: realization

  write(realization%option%io_buffer, &
        '(/,"PFLOTRAN must be compiled with -DUSE_HDF5 to ", &
        &"read HDF5 formatted structured grids.",/)')
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

  interface

     subroutine SAMRCopyVecToVecComponent(vec,svec, component)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       Vec :: vec, svec
       PetscInt :: component
     end subroutine SAMRCopyVecToVecComponent

     subroutine SAMRRegisterForViz(ptr,vec,component,dname,dnamec)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: ptr
       Vec :: vec
       PetscInt :: component
       PetscInt :: dname, dnamec
     end subroutine SAMRRegisterForViz

     subroutine SAMRWritePlotData(ptr, time)
#include "finclude/petsc.h"
       PetscFortranAddr :: ptr
       PetscReal :: time
     end subroutine SAMRWritePlotData

  end interface

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
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscTruth :: first

  discretization => realization%discretization
  patch => realization%patch
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  output_option => realization%output_option

  first = hdf5_first
  if (option%restart_flag .and. option%restart_time > 1.d-40)  &
    first = PETSC_FALSE

  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.h5'

  if(.not.(option%use_samr)) then
     
     grid => patch%grid
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
        string = "Coordinates"
        call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

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

        call h5gclose_f(grp_id,hdf5_err)

     endif

     ! create a group for the data set
     write(string,'(''Time:'',es12.4,x,a1)') &
          option%time/output_option%tconv,output_option%tunit
     if (len_trim(output_option%plot_name) > 2) then
        string = trim(string) // ' ' // output_option%plot_name
     endif
     call h5eset_auto_f(OFF,hdf5_err)
     call h5gopen_f(file_id,string,grp_id,hdf5_err)
     if (hdf5_err /= 0) then
       call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
     endif
     call h5eset_auto_f(ON,hdf5_err)
  
  else

     app_ptr = discretization%amrgrid%p_application

     select case(option%iflowmode)
        case(RICHARDS_MODE)
           nviz_flow = 2
        case(MPH_MODE)
           nviz_flow = 7+2*option%nflowspec
        case(IMS_MODE)
           nviz_flow = 2+4*option%nphase
        case(THC_MODE)
           nviz_flow = 4+option%nflowspec
        case default
           nviz_flow=0
     end select

     if (option%ntrandof > 0) then
        if (associated(reaction)) then
           nviz_tran = option%ntrandof+reaction%nkinmnrl
        else
           nviz_tran = reaction%nkinmnrl
        endif
     endif

     nviz_dof = nviz_flow+nviz_tran

!!$     if (associated(patch%imat)) then
!!$        nviz_dof = nviz_dof+1
!!$     endif

     if(first) then
        call AMRGridCreateVector(discretization%amrgrid, nviz_dof, field%samr_viz_vec, &
             GLOBAL, PETSC_TRUE, option)
     endif

     current_component = 0
  endif

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)

  select case(option%iflowmode)
  
    case(MPH_MODE,THC_MODE, IMS_MODE,&
         RICHARDS_MODE)

      ! temperature
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Temperature"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,TEMPERATURE,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select

      ! pressure
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Pressure"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,PRESSURE,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select

      ! liquid saturation
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,RICHARDS_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Liquid Saturation"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)  
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,LIQUID_SATURATION,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select

      ! gas saturation
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Gas Saturation"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,GAS_SATURATION,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select
      
      ! liquid energy
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Liquid Energy"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,LIQUID_ENERGY,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select
      
      ! gas energy
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE)    
          call OutputGetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          if(.not.(option%use_samr)) then
             string = "Gas Energy"
             call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
             call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,GAS_ENERGY,ZERO_INTEGER)
             endif
             current_component=current_component+1
          endif
      end select
    
      ! liquid mole fractions
      select case(option%iflowmode)
        case (MPH_MODE,THC_MODE,IMS_MODE)
          do i=1,option%nflowspec
            call OutputGetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i)
            if(.not.(option%use_samr)) then
               write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
               call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
            else
               call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,LIQUID_MOLE_FRACTION,i)
             endif
               current_component=current_component+1
            endif
          enddo
      end select
      
      ! gas mole fractions
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE)      
          do i=1,option%nflowspec
             call OutputGetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i)
             if(.not.(option%use_samr)) then
                write(string,'(''Gas Mole Fraction('',i4,'')'')') i
                call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
            else
               call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
             if(first) then
                call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,GAS_MOLE_FRACTION,i)
             endif
               current_component=current_component+1
            endif
          enddo
      end select
#if 0    
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,ZERO_INTEGER)
        string = "Volume Fraction"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      endif
#endif    
      ! phase
      select case(option%iflowmode)
        case (MPH_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          if (.not.(option%use_samr)) then
            string = "Phase"
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,PHASE,ZERO_INTEGER)
            endif
            current_component=current_component+1
          endif
                
      end select
  
    case default

  end select

  if (option%ntrandof > 0) then
    if (associated(reaction)) then
      if (reaction%print_pH .and. reaction%h_ion_id > 0) then
        call OutputGetVarFromArray(realization,global_vec,PH,reaction%h_ion_id)
        if (.not.(option%use_samr)) then
          write(string,'(''pH'')')
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
        else
          call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
          if(first) then
             call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,PH,reaction%h_ion_id)
          endif
          current_component=current_component+1
        endif
      endif
      do i=1,reaction%ncomp
        if (reaction%primary_species_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,TOTAL_MOLARITY,i)
          if (.not.(option%use_samr)) then
            write(string,'(a)') reaction%primary_species_names(i)
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,TOTAL_MOLARITY,i)
            endif
            current_component=current_component+1
          endif
        endif
      enddo
      if (realization%output_option%print_act_coefs) then
        do i=1,reaction%ncomp
          if (reaction%primary_species_print(i)) then
            call OutputGetVarFromArray(realization,global_vec,PRIMARY_ACTIVITY_COEF,i)
            if (.not.(option%use_samr)) then
              write(string,'(a)') trim(reaction%primary_species_names(i)) // '_gam'
              call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
            else
              call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
              if(first) then
                 call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,PRIMARY_ACTIVITY_COEF,i)
              endif
              current_component=current_component+1
            endif
          endif
        enddo
      endif
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_VOLUME_FRACTION,i)
          if (.not.(option%use_samr)) then
            write(string,'(a)') trim(reaction%kinmnrl_names(i)) // '_vf'
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,MINERAL_VOLUME_FRACTION,i)
            endif
            current_component=current_component+1
          endif
        endif
      enddo
      do i=1,reaction%nkinmnrl
        if (reaction%kinmnrl_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,MINERAL_RATE,i)
          if (.not.(option%use_samr)) then
            write(string,'(a)') trim(reaction%kinmnrl_names(i)) // '_rt'
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,MINERAL_RATE,i)
            endif
            current_component=current_component+1
          endif
        endif
      enddo
      do i=1,reaction%neqsurfcmplxrxn
        if (reaction%surface_site_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX_FREE,i)
          if (.not.(option%use_samr)) then
            write(string,'(a)') reaction%surface_site_names(i)
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,SURFACE_CMPLX_FREE,i)
            endif
            current_component=current_component+1
          endif
        endif
      enddo
      do i=1,reaction%neqsurfcmplx
        if (reaction%surface_complex_print(i)) then
          call OutputGetVarFromArray(realization,global_vec,SURFACE_CMPLX,i)
          if (.not.(option%use_samr)) then
            write(string,'(a)') reaction%surface_complex_names(i)
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
          else
            call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
            if(first) then
               call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,SURFACE_CMPLX,i)
            endif
            current_component=current_component+1
          endif
        endif
      enddo
    endif
  endif
  
  ! material id

  if(.not.(option%use_samr)) then
     if(associated(patch%imat)) then
        call OutputGetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
        if(.not.(option%use_samr)) then
           string = "Material_ID"
           call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
        else
           call SAMRCopyVecToVecComponent(global_vec,field%samr_viz_vec, current_component)
           if(first) then
              call SAMRRegisterForViz(app_ptr,field%samr_viz_vec,current_component,MATERIAL_ID,ZERO_INTEGER)
       endif
       current_component=current_component+1
    endif
    endif
  endif
  
  if(.not.(option%use_samr)) then
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
  endif  
  ! call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  if(.not.(option%use_samr)) then
     call h5gclose_f(grp_id,hdf5_err)
     call h5fclose_f(file_id,hdf5_err)
     call h5close_f(hdf5_err)
  else
     call SAMRWritePlotData(app_ptr, option%time/output_option%tconv)
!     call VecDestroy(field%samr_viz_vec,ierr)
  endif
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
  
  use AMR_Grid_Module
 
#ifndef USE_HDF5
  implicit none
  
  type(realization_type) :: realization

  write(realization%option%io_buffer, &
        '(/,"PFLOTRAN must be compiled with -DUSE_HDF5 to ", &
        &"read HDF5 formatted structured grids.",/)')
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
  Vec :: samr_vec
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
  call OutputGetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
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

  PetscTruth, save :: trick_flux_vel_x = PETSC_FALSE
  PetscTruth, save :: trick_flux_vel_y = PETSC_FALSE
  PetscTruth, save :: trick_flux_vel_z = PETSC_FALSE

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
    if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER,MPI_INTEGER,MPI_MIN,option%mycomm,ierr)
    if (i == 0) trick_flux_vel_x = PETSC_TRUE
    if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER,MPI_INTEGER,MPI_MIN,option%mycomm,ierr)
    if (j == 0) trick_flux_vel_y = PETSC_TRUE
    if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER,MPI_INTEGER,MPI_MIN,option%mycomm,ierr)
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
      if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = PETSC_TRUE
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = PETSC_TRUE
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
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

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                        nx_global,ny_global,nz_global, &
                        nx_local,ny_local,nz_local, &
                        grid%structured_grid%nxs,grid%structured_grid%nys,grid%structured_grid%nzs)
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
  if (option%myrank == option%io_rank) then
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
                                      num_cells,ghosted_ids)

  use Realization_module
  use Grid_module
  use Option_module

  implicit none
  
  PetscReal :: OutputGetVarFromArrayAtCoord
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
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
    value = RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id)
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
subroutine OutputGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  call PetscLogEventBegin(logging%event_output_get_var_from_array, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
                        
  call RealizationGetDataset(realization,vec,ivar,isubvar)

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
                          patch%boundary_velocities(1,sum_connection)* &
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
! OutputMassBalance: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 06/18/08
!
! ************************************************************************** !  
#if 0
subroutine OutputMassBalance(realization)

  use Realization_module
  use Discretization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i, comma_count, quote_count
  character(len=MAXSTRINGLENGTH) :: filename
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
  else
    if (output_option%plot_number < 10) then
      write(string,'("00",i1,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(string,'("0",i2,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(string,'(i3,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(string,'(i4,".tec")') output_option%plot_number  
    endif
    filename = 'mass-' // trim(option%group_prefix) // '-' //trim(string)
  endif
  
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot output file: ' // trim(filename)
    call printMsg(option)    
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
      case (IMS_MODE)
      case(THC_MODE)
        string2 = '"Water [kg]","Solute [mol]","Energy []"'
      case(RICHARDS_MODE)
    end select
    string = trim(string) // trim(string2)
    
    ! write transport variables
    if (option%ntrandof > 0) then
 !     string2 = RTGetTecplotHeader(realization)
 !     string = trim(string) // trim(string2)
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
    case(THC_MODE)
#if 0
      ! water mass
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! solute mass
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! energy
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,IMS_MODE)
          call OutputGetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
#endif
  end select

end subroutine OutputMassBalance
#endif

! ************************************************************************** !
!
! OutputMassBalanceNew: Print to Tecplot POINT format
! author: Glenn Hammond
! date: 06/18/08
!
! ************************************************************************** !  
subroutine OutputMassBalanceNew(realization)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  
  use Richards_module
  use Reactive_Transport_module
  
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
  type(coupler_type), pointer :: boundary_condition
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: word
  character(len=4) :: strcol
  PetscInt :: fid = 86
  PetscInt :: ios
  PetscInt :: i,icol
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: offset
  PetscInt :: iphase
  PetscInt :: icomp
  PetscReal :: sum_kg(realization%option%nphase)
  PetscReal :: sum_kg_global(realization%option%nphase)
  PetscReal :: sum_mol(realization%option%ntrandof,realization%option%nphase)
  PetscReal :: sum_mol_global(realization%option%ntrandof,realization%option%nphase)
  PetscTruth :: local_first
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  reaction => realization%reaction
  output_option => realization%output_option
  
  local_first = mass_balance_first
  
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.dat'
  else
    filename = 'mass_balance' // trim(option%group_prefix) // '.dat'
  endif
  
  ! open file
  if (option%myrank == option%io_rank) then

    option%io_buffer = '--> write tecplot mass balance file: ' // trim(filename)
    call printMsg(option)    

    if (mass_balance_first .and. option%restart_flag) then ! check if file already exists
      ios = 0
      open(unit=fid,file=filename,action="write",status="old",iostat=ios)
      if (ios == 0) then
        close(fid)
        local_first = PETSC_FALSE 
      endif
    endif

    if (local_first) then
      open(unit=fid,file=filename,action="write",status="replace")

      ! write header
      write(fid,'(a)',advance="no") '"Time[' // trim(output_option%tunit) // ']"'  
      
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          write(fid,'(a)',advance="no") ',"2-Global Water Mass [kg]"'
      end select
      
      icol = 2
      if (option%ntrandof > 0) then
        do i=1,reaction%ncomp
          if (reaction%primary_species_print(i)) then
            icol = icol + 1
            write(strcol,'(i3,"-")') icol
            write(fid,'(a)',advance="no") ',"' // trim(strcol) // 'Global ' // &
                trim(reaction%primary_species_names(i)) // ' [mol]"'
          endif
        enddo
      endif
      
      boundary_condition => patch%boundary_conditions%first
      do 
        if (.not.associated(boundary_condition)) exit

        select case(option%iflowmode)
          case(RICHARDS_MODE)
            icol = icol + 1
            write(strcol,'(i3,"-")') icol
            write(fid,'(a)',advance="no") ',"' // trim(strcol) // &
              trim(boundary_condition%name) // ' Water Mass [kg]"'
            icol = icol + 1
            write(strcol,'(i3,"-")') icol
            write(fid,'(a)',advance="no") ',"' // trim(strcol) // &
              trim(boundary_condition%name) // ' Water Mass [kg/' // &
              trim(output_option%tunit) // ']"'
        end select
        
        if (option%ntrandof > 0) then
          do i=1,reaction%ncomp
            if (reaction%primary_species_print(i)) then
              icol = icol + 1
              write(strcol,'(i3,"-")') icol
              write(fid,'(a)',advance="no") ',"' // trim(strcol) // &
                  trim(boundary_condition%name) // ' ' // &
                  trim(reaction%primary_species_names(i)) // ' [mol]"'
            endif
          enddo

          do i=1,reaction%ncomp
            if (reaction%primary_species_print(i)) then
              icol = icol + 1
              write(strcol,'(i3,"-")') icol
              write(fid,'(a)',advance="no") ',"' // trim(strcol) // &
                  trim(boundary_condition%name) // ' ' // &
                  trim(reaction%primary_species_names(i)) // ' [mol/' // &
                  trim(output_option%tunit) // ']"'
            endif
          enddo
        endif
        boundary_condition => boundary_condition%next
      
      enddo
      
#ifdef COMPUTE_INTERNAL_MASS_FLUX
      do offset = 1, 4
        write(word,'(i6)') offset*100
        select case(option%iflowmode)
          case(RICHARDS_MODE)
            write(fid,'(a)',advance="no") ',"' // &
              trim(adjustl(word)) // 'm Water Mass [kg]"'
        end select
        
        if (option%ntrandof > 0) then
          do i=1,reaction%ncomp
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

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  if (option%ntrandof > 0) then
    rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  endif    
  
  if (option%nflowdof > 0) then
    sum_kg = 0.d0
    call RichardsComputeMassBalance(realization,sum_kg)
    call MPI_Reduce(sum_kg,sum_kg_global, &
                    option%nphase,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)
                        
    if (option%myrank == option%io_rank) then
      write(fid,110,advance="no") sum_kg_global
    endif
  endif
  
  if (option%ntrandof > 0) then
    sum_mol = 0.d0
    call RTComputeMassBalance(realization,sum_mol)
    call MPI_Reduce(sum_mol,sum_mol_global,option%nphase*reaction%ncomp, &
                    MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

    if (option%myrank == option%io_rank) then
      do iphase = 1, option%nphase
        do icomp = 1, reaction%ncomp
          if (reaction%primary_species_print(icomp)) then
            write(fid,110,advance="no") sum_mol_global(icomp,iphase)
          endif
        enddo
      enddo
    endif
  endif
  
  boundary_condition => patch%boundary_conditions%first
  do 
    if (.not.associated(boundary_condition)) exit

    offset = boundary_condition%connection_set%offset
    
    if (option%nflowdof > 0) then
      ! print out cumulative H2O flux
      sum_kg = 0.d0
      do iconn = 1, boundary_condition%connection_set%num_connections
        sum_kg = sum_kg + global_aux_vars_bc(offset+iconn)%mass_balance
      enddo

      call MPI_Reduce(sum_kg,sum_kg_global, &
                      option%nphase,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                          
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        write(fid,110,advance="no") -sum_kg_global
      endif

      ! print out H2O flux
      sum_kg = 0.d0
      do iconn = 1, boundary_condition%connection_set%num_connections
        sum_kg = sum_kg + global_aux_vars_bc(offset+iconn)%mass_balance_delta
      enddo

      call MPI_Reduce(sum_kg,sum_kg_global, &
                      option%nphase,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                          
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        write(fid,110,advance="no") -sum_kg_global*output_option%tconv
      endif
    endif
    
    if (option%ntrandof > 0) then

      ! print out cumulative boundary flux
      sum_mol = 0.d0
      do iconn = 1, boundary_condition%connection_set%num_connections
        sum_mol = sum_mol + rt_aux_vars_bc(offset+iconn)%mass_balance
      enddo

      call MPI_Reduce(sum_mol,sum_mol_global,option%nphase*reaction%ncomp, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        do iphase = 1, option%nphase
          do icomp = 1, reaction%ncomp
            if (reaction%primary_species_print(icomp)) then
              write(fid,110,advance="no") -sum_mol_global(icomp,iphase)
            endif
          enddo
        enddo
      endif
    
      ! print out boundary flux
      sum_mol = 0.d0
      do iconn = 1, boundary_condition%connection_set%num_connections
        sum_mol = sum_mol + rt_aux_vars_bc(offset+iconn)%mass_balance_delta 
      enddo

      call MPI_Reduce(sum_mol,sum_mol_global,option%nphase*reaction%ncomp, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                      
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        do iphase = 1, option%nphase
          do icomp = 1, reaction%ncomp
            if (reaction%primary_species_print(icomp)) then
              write(fid,110,advance="no") -sum_mol_global(icomp,iphase)* &
                                          output_option%tconv
            endif
          enddo
        enddo
      endif
    endif

    boundary_condition => boundary_condition%next
  
  enddo

#ifdef COMPUTE_INTERNAL_MASS_FLUX

  do offset = 1, 4
    iconn = offset*20-1

    if (option%nflowdof > 0) then
      sum_kg = 0.d0
      sum_kg = sum_kg + patch%aux%Global%aux_vars(iconn)%mass_balance

      call MPI_Reduce(sum_kg,sum_kg_global, &
                      option%nphase,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)
                          
      if (option%myrank == option%io_rank) then
        ! change sign for positive in / negative out
        write(fid,110,advance="no") -sum_kg_global
      endif
    endif
    
    if (option%ntrandof > 0) then

      sum_mol = 0.d0
      sum_mol = sum_mol + patch%aux%RT%aux_vars(iconn)%mass_balance

      call MPI_Reduce(sum_mol,sum_mol_global,option%nphase*reaction%ncomp, &
                      MPI_DOUBLE_PRECISION,MPI_SUM, &
                      option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        do iphase = 1, option%nphase
          do icomp = 1, reaction%ncomp
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

end subroutine OutputMassBalanceNew

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

end module Output_module

    
       
