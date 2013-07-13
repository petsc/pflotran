#ifdef SURFACE_FLOW

module Surface_Init_module

  implicit none

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"


  public :: SurfaceInitReadRequiredCards, &
            SurfaceInitReadInput, &
            SurfaceInitMatPropToRegions, &
            SurfaceInitReadRegionFiles
contains

! ************************************************************************** !
!> This routine reads the required input file cards related to surface flows
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/18/12
! ************************************************************************** !

subroutine SurfaceInitReadRequiredCards(surf_realization)

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Level_module

  use Surface_Realization_class
  use Surface_Auxiliary_module

  implicit none

  type(surface_realization_type)     :: surf_realization
  type(discretization_type), pointer :: discretization

  character(len=MAXSTRINGLENGTH) :: string
  
  type(patch_type), pointer   :: patch
  type(level_type), pointer   :: level
  type(grid_type), pointer    :: grid
  type(option_type), pointer  :: option
  type(input_type), pointer   :: input
  
  patch          => surf_realization%patch
  option         => surf_realization%option
  discretization => surf_realization%discretization
  
  input => surf_realization%input
  
! Read in select required cards
!.........................................................................
 
  ! GRID information
!  string = "GRID"
!  call InputFindStringInFile(input,option,string)
!  call InputFindStringErrorMsg(input,option,string)

  ! SURFACE_FLOW information
  string = "SURFACE_FLOW"
  call InputFindStringInFile(input,option,string)
  if(InputError(input)) return
  option%nsurfflowdof = 1
  
  string = "SURF_GRID"
  call InputFindStringInFile(input,option,string)
!  call SurfaceFlowReadRequiredCardsFromInput(surf_realization,input,option)
  call SurfaceInit(surf_realization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      patch%surf_or_subsurf_flag = SURFACE
      if (.not.associated(surf_realization%level_list)) then
        surf_realization%level_list => LevelCreateList()
      endif
      level => LevelCreate()
      call LevelAddToList(level,surf_realization%level_list)
      call PatchAddToList(patch,level%patch_list)
      surf_realization%patch => patch
  end select
    
end subroutine SurfaceInitReadRequiredCards

! ************************************************************************** !
!> This routine reads required surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceInit(surf_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Surface_Material_module
  use Surface_Realization_class
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Discretization_module
  use Region_module
  use Condition_module
  use Unstructured_Grid_Aux_module

  implicit none

  type(surface_realization_type)               :: surf_realization
  type(discretization_type),pointer            :: discretization
  type(grid_type), pointer                     :: grid
  type(input_type)                             :: input
  type(option_type)                            :: option
  type(unstructured_grid_type), pointer        :: un_str_sfgrid
  character(len=MAXWORDLENGTH)                 :: word
  character(len=MAXWORDLENGTH)                 :: unstructured_grid_ctype
  PetscInt :: unstructured_grid_itype

  discretization => surf_realization%discretization

  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  call InputReadFlotranString(input,option)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
  call StringToUpper(word)
    
  select case(trim(word))
    case ('TYPE')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword','TYPE')
      call StringToUpper(word)

      select case(trim(word))
        case ('UNSTRUCTURED')
          unstructured_grid_itype = IMPLICIT_UNSTRUCTURED_GRID
          unstructured_grid_ctype = 'implicit unstructured'
          discretization%itype = UNSTRUCTURED_GRID
          call InputReadNChars(input,option, &
                               discretization%filename, &
                               MAXSTRINGLENGTH, &
                               PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','filename')

          grid => GridCreate()
          un_str_sfgrid => UGridCreate()
          un_str_sfgrid%grid_type = TWO_DIM_GRID
          if (index(discretization%filename,'.h5') > 0) then
            call UGridReadHDF5SurfGrid( un_str_sfgrid, &
                                        discretization%filename, &
                                        option)
          else
            call UGridReadSurfGrid(un_str_sfgrid, &
                                   surf_realization%subsurf_filename, &
                                   discretization%filename, &
                                   option)
          endif
          grid%unstructured_grid => un_str_sfgrid
          discretization%grid => grid
          grid%itype = unstructured_grid_itype
          grid%ctype = unstructured_grid_ctype

        case default
          option%io_buffer = 'Surface-flow supports only unstructured grid'
          call printErrMsg(option)
      end select
  end select

end subroutine SurfaceInit

! ************************************************************************** !
!> This routine reads surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceInitReadInput(surf_realization,surf_flow_solver,input,option)

  use Option_module
  use Input_module
  use String_module
  use Surface_Material_module
  use Surface_Realization_class
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Unstructured_Grid_Aux_module
  use Discretization_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module
  use Debug_module
  use Units_module
  use Waypoint_module
  use Patch_module
  use Solver_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Output_Surface_module

  implicit none

  type(surface_realization_type)               :: surf_realization
  type(solver_type)                            :: surf_flow_solver
  type(input_type)                             :: input
  type(option_type)                            :: option
  
  type(discretization_type),pointer            :: discretization
  type(grid_type), pointer                     :: grid
  type(unstructured_grid_type), pointer        :: un_str_sfgrid
  type(surface_material_property_type),pointer :: surf_material_property
  type(region_type), pointer                   :: region
  type(flow_condition_type), pointer           :: flow_condition
  type(coupler_type), pointer                  :: coupler
  type(strata_type), pointer                   :: strata
  class(dataset_base_type), pointer            :: dataset

  type(patch_type), pointer                    :: patch
  type(output_option_type), pointer            :: output_option
  PetscReal :: units_conversion

  character(len=MAXWORDLENGTH)                 :: word
  character(len=MAXWORDLENGTH)                 :: card
  character(len=1) :: backslash

  PetscBool :: velocities
  PetscBool :: fluxes
  PetscBool :: continuation_flag
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate

  type(waypoint_type), pointer :: waypoint
  PetscReal :: temp_real, temp_real2

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++

  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  discretization => surf_realization%discretization
  output_option => surf_realization%output_option

  patch => surf_realization%patch

  if (associated(patch)) grid => patch%grid

  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
    call StringToUpper(word)
    write(*,*),'word :: ',trim(word)

    select case(trim(word))
      !.........................................................................
      ! Read surface grid information
      case ('SURF_GRID')
        call InputSkipToEND(input,option,trim(word))
      !.........................................................................
      case ('SURF_FLOW_FORMULATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call StringToUpper(word)
        select case(trim(word))
          case ('KINEMATIC')
            option%surface_flow_formulation = KINEMATIC_WAVE
          case ('DIFFUSIVE')
            option%surface_flow_formulation = DIFFUSION_WAVE
          case default
            option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
              'not recognized'
            call printErrMsg(option)
        end select

      !.........................................................................
      ! Read surface material information
      case ('SURF_MATERIAL_PROPERTY')
        surf_material_property => SurfaceMaterialPropertyCreate()

        call InputReadWord(input,option,surf_material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        call SurfaceMaterialPropertyRead(surf_material_property,input,option)
        call SurfaceMaterialPropertyAddToList(surf_material_property, &
                                          surf_realization%surf_material_properties)
        nullify(surf_material_property)

      !.........................................................................
      case ('SURF_REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SURF_REGION')
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,surf_realization%surf_regions)
        nullify(region)

      !.........................................................................
      case ('SURF_FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_FLOW_CONDITION','name')
        call printMsg(option,flow_condition%name)
        if (option%iflowmode == G_MODE) then
          call FlowConditionGeneralRead(flow_condition,input,option)
        else
          call FlowConditionRead(flow_condition,input,option)
        endif
        call FlowConditionAddToList(flow_condition,surf_realization%surf_flow_conditions)
        nullify(flow_condition)

      !.........................................................................
      case ('SURF_BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name')
        call CouplerRead(coupler,input,option)
        call SurfRealizAddCoupler(surf_realization,coupler)
        nullify(coupler)

      !.........................................................................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call SurfRealizAddStrata(surf_realization,strata)
        nullify(strata)
        
      !.........................................................................
      case ('SURF_INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call SurfRealizAddCoupler(surf_realization,coupler)
        nullify(coupler)        

      !.........................................................................
      case ('SURF_SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call CouplerRead(coupler,input,option)
        call SurfRealizAddCoupler(surf_realization,coupler)
        nullify(coupler)        

      !.........................................................................
      case ('SURF_SUBSURFACE_COUPLING')
        call InputReadFlotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call StringToUpper(word)
        select case(trim(word))
          case('DECOUPLED')
            option%subsurf_surf_coupling = DECOUPLED
          case('SEQ_COUPLED')
            option%subsurf_surf_coupling = SEQ_COUPLED
          case('FULLY_COUPLED')
            option%subsurf_surf_coupling = FULLY_COUPLED
          case default
            option%io_buffer = 'Invalid value for SURF_SUBSURFACE_COUPLING'
            call printErrMsg(option)
        end select
        call InputSkipToEND(input,option,trim(word))

      !.........................................................................
      case ('SURF_DEBUG')
        call DebugRead(surf_realization%debug,input,option)

      !.........................................................................
      case ('SURF_OUTPUT')
        velocities = PETSC_FALSE
        fluxes = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','SURF_OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial = PETSC_FALSE
            case('PERMEABILITY')
              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              output_option%print_porosity = PETSC_TRUE
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE
            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','SURF_OUTPUT')
              units_conversion = UnitsConvertToInternal(word,option)
              continuation_flag = PETSC_TRUE
              do
                continuation_flag = PETSC_FALSE
                if (index(input%buf,backslash) > 0) &
                  continuation_flag = PETSC_TRUE
                input%ierr = 0
                do
                  if (InputError(input)) exit
                  call InputReadDouble(input,option,temp_real)
                  if (.not.InputError(input)) then
                    waypoint => WaypointCreate()
                    waypoint%time = temp_real*units_conversion
                    waypoint%print_output = PETSC_TRUE
                    write(*,*),'Inserting waypoint in surf_realization: ',waypoint%time
                    call WaypointInsertInList(waypoint,surf_realization%waypoints)
                  endif
                enddo
                if (.not.continuation_flag) exit
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
              enddo
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'SURF_OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,OUTPUT_FILE')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,OUTPUT_FILE.'
                  call printErrMsg(option)
              end select
            case('SCREEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_screen = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,SCREEN')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,SCREEN.'
                  call printErrMsg(option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'SURF_OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'SURF_OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'SURF_OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word,option)
                  output_option%periodic_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then

                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word,option)
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'SURF_OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_output = PETSC_TRUE
                        write(*,*),'Inserting waypoint in surf_realization: >>>>>>>> ',waypoint%time
                        call WaypointInsertInList(waypoint,surf_realization%waypoints)
                        temp_real = temp_real + output_option%periodic_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'SURF_OUTPUT,PERIODIC,TIME')
                    endif
                  endif
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,PERIODIC,'// &
                                     'TIMESTEP.'
                  call printErrMsg(option)
              end select
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                'OUTPUT, PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  units_conversion = UnitsConvertToInternal(word,option) 
                  output_option%periodic_tr_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_tr_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,'// &
                                     'PERIODIC_OBSERVATION,TIMESTEP.'
                  call printErrMsg(option)
              end select
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 1) then 
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recongnized.  Use "SINGLE_FILE" or ' // &
                          '"MULTIPLE_FILES".'
                        call printErrMsg(option)
                    end select
                  endif
                case ('MAD')
                  output_option%print_mad = PETSC_TRUE
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'TECPLOT','OUTPUT,FORMAT') 
                  call StringToUpper(word)
                  select case(trim(word))
                    case('POINT')
                      output_option%tecplot_format = TECPLOT_POINT_FORMAT
                    case('BLOCK')
                      output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                    case('FEQUADRILATERAL')
                      output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT
                    case default
                      option%io_buffer = 'TECPLOT format (' // trim(word) // &
                                         ') not recongnized.'
                      call printErrMsg(option)
                  end select
                  if (output_option%tecplot_format == TECPLOT_POINT_FORMAT &
                      .and. option%mycommsize > 1) then
                    output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                  endif
                  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
                    output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,FORMAT.'
                  call printErrMsg(option)
              end select

            case('VELOCITIES')
              velocities = PETSC_TRUE
            case('FLUXES')
              fluxes = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
            case('HYDROGRAPH')
              output_option%print_hydrograph = PETSC_TRUE
            case('PROCESSOR_ID')
              output_option%print_iproc = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case('VARIABLES')
              call OutputSurfaceVariableRead(input,option,output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputSurfaceVariableRead(input,option,output_option%aveg_output_variable_list)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in OUTPUT.'
              call printErrMsg(option)
          end select
        enddo

        if (velocities) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_velocities = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_velocities = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_velocities = PETSC_TRUE
        endif
        if (fluxes) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_flux_velocities = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_flux_velocities = PETSC_TRUE
        endif
        if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate.or.aveg_energy_flowrate) then
          if (output_option%print_hdf5) then
#ifndef STORE_FLOWRATES
            option%io_buffer='To output FLOWRATES/MASS_FLOWRATE/ENERGY_FLOWRATE, '// &
              'compile with -DSTORE_FLOWRATES'
            call printErrMsg(option)
#endif
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if(aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if(output_option%periodic_output_time_incr==0.d0) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ ' // &
                  'AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without' // &
                  ' PERIODIC TIME being set.'
                call printErrMsg(option)
              endif
            endif
           option%store_flowrate = PETSC_TRUE
          else
            option%io_buffer='Output FLOWRATES/MASS_FLOWRATE/ENERGY_FLOWRATE ' // &
              'only available in HDF5 format'
            call printErrMsg(option)
          endif
        endif

      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadNewton(surf_flow_solver,input,option)
        end select

      !.........................................................................
      case ('SURF_TIME')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','SURF_TIME')
          select case(trim(word))
            case ('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Initial Timestep Size Time Units','TIME')
              surf_realization%dt_min = temp_real*UnitsConvertToInternal(word,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Maximum Timestep Size Time Units','TIME')
              surf_realization%dt_max = temp_real*UnitsConvertToInternal(word,option)
            case('COUPLING_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Coupling Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Coupling Timestep Size Time Units','TIME')
              surf_realization%dt_coupling = temp_real*UnitsConvertToInternal(word,option)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in TIME.'
              call printErrMsg(option)
            end select
        enddo
      !.........................................................................
      case ('SURF_DATASET')
      nullify(dataset)
      call DatasetRead(input,dataset,option)
      call DatasetBaseAddToList(dataset,surf_realization%datasets)
      nullify(dataset)
      
      !.........................................................................
      case ('SURF_RESTART')
        option%surf_restart_flag = PETSC_TRUE
        call InputReadNChars(input,option,option%surf_restart_filename,MAXSTRINGLENGTH, &
                             PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_RESTART','Surfae restart file name') 
        call InputReadDouble(input,option,option%restart_time)
        if (input%ierr == 0) then
          call printErrMsg(option,'Setting time to value not supported in surface-flow')
        endif

      case default
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select
  enddo

  if(option%restart_flag .neqv. option%surf_restart_flag) then
    option%io_buffer='option%restart_flag /= option%surf_restart_flag'
    call printErrMsg(option)
  endif

end subroutine SurfaceInitReadInput

! ************************************************************************** !
!> This routine assigns surface material properties to associated regions in
!! the model (similar to assignMaterialPropToRegions)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/13/12
! ************************************************************************** !

subroutine SurfaceInitMatPropToRegions(surf_realization)

  use Surface_Realization_class
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Level_module
  use Surface_Field_module
  use Surface_Material_module
  
  use HDF5_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  PetscReal, pointer :: man0_p(:)
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, surf_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: patch  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  type(surface_material_property_type), pointer :: surf_material_property
  type(surface_material_property_type), pointer :: null_surf_material_property
  type(region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  
  option => surf_realization%option
  discretization => surf_realization%discretization
  surf_field => surf_realization%surf_field

  ! loop over all patches and allocation material id arrays
  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      if (.not.associated(cur_patch%imat)) then
        allocate(cur_patch%imat(cur_patch%grid%ngmax))
        ! initialize to "unset"
        cur_patch%imat = -999
        ! also allocate saturation function id
        allocate(cur_patch%sat_func_id(cur_patch%grid%ngmax))
        cur_patch%sat_func_id = -999
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      grid => cur_patch%grid
      strata => cur_patch%strata%first
      do
        if (.not.associated(strata)) exit
        ! Read in cell by cell material ids if they exist
        if (.not.associated(strata%region) .and. strata%active) then
          option%io_buffer = 'Reading of material prop from file for' // &
            ' surface flow is not implemented.'
          call printErrMsgByRank(option)
          !call readMaterialsFromFile(realization,strata%realization_dependent, &
          !                           strata%material_property_filename)
        ! Otherwise, set based on region
        else if (strata%active) then
          update_ghosted_material_ids = PETSC_TRUE
          region => strata%region
          surf_material_property => strata%surf_material_property
          if (associated(region)) then
            istart = 1
            iend = region%num_cells
          else
            istart = 1
            iend = grid%nlmax
          endif
          do icell=istart, iend
            if (associated(region)) then
              local_id = region%cell_ids(icell)
            else
              local_id = icell
            endif
            ghosted_id = grid%nL2G(local_id)
            cur_patch%imat(ghosted_id) = surf_material_property%id
          enddo
        endif
        strata => strata%next
      enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call SurfRealizLocalToLocalWithArray(surf_realization,MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_surf_material_property => SurfaceMaterialPropertyCreate()
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      call GridVecGetArrayF90(grid,surf_field%mannings0,man0_p,ierr)

      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        surf_material_id = cur_patch%imat(ghosted_id)
        if (surf_material_id == 0) then ! accomodate inactive cells
          surf_material_property = null_surf_material_property
        else if ( surf_material_id > 0 .and. &
                  surf_material_id <= &
                  size(surf_realization%surf_material_property_array)) then
          surf_material_property => &
            surf_realization%surf_material_property_array(surf_material_id)%ptr
          if (.not.associated(surf_material_property)) then
            write(dataset_name,*) surf_material_id
            option%io_buffer = 'No material property for surface material id ' // &
                               trim(adjustl(dataset_name)) &
                               //  ' defined in input file.'
            call printErrMsgByRank(option)
          endif
        else if (surf_material_id < -998) then 
          write(dataset_name,*) grid%nG2A(ghosted_id)
          option%io_buffer = 'Uninitialized surface material id in patch at cell ' // &
                             trim(adjustl(dataset_name))
          call printErrMsgByRank(option)
        else if (surf_material_id > size(surf_realization%surf_material_property_array)) then
          write(option%io_buffer,*) surf_material_id
          option%io_buffer = 'Unmatched surface material id in patch:' // &
            adjustl(trim(option%io_buffer))
          call printErrMsgByRank(option)
        else
          option%io_buffer = 'Something messed up with surface material ids. ' // &
            ' Possibly material ids not assigned to all grid cells. ' // &
            ' Contact Glenn!'
          call printErrMsgByRank(option)
        endif
        man0_p(local_id) = surf_material_property%mannings
      enddo ! local_id - loop

      call GridVecRestoreArrayF90(grid,surf_field%mannings0,man0_p,ierr)
      
      cur_patch => cur_patch%next
    enddo ! looping over patches
    cur_level => cur_level%next
  enddo ! looping over levels
  
  call SurfaceMaterialPropertyDestroy(null_surf_material_property)
  nullify(null_surf_material_property)

  call DiscretizationGlobalToLocal(discretization,surf_field%mannings0, &
                                   surf_field%mannings_loc,ONEDOF)

end subroutine SurfaceInitMatPropToRegions

! ************************************************************************** !
!> This routine reads surface region files
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/20/12
! ************************************************************************** !
subroutine SurfaceInitReadRegionFiles(surf_realization)

  use Surface_Realization_class
  use Region_module
  use HDF5_module
  use Grid_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(region_type), pointer :: surf_region
  
  surf_region => surf_realization%surf_regions%first
  do 
    if (.not.associated(surf_region)) exit
    if (len_trim(surf_region%filename) > 1) then
      if (index(surf_region%filename,'.h5') > 0) then
        if (surf_region%grid_type == STRUCTURED_GRID) then
          !call HDF5ReadRegionFromFile(surf_realization,surf_region,surf_region%filename)
        else
#if defined(PETSC_HAVE_HDF5)
          call HDF5ReadUnstructuredGridRegionFromFile(surf_realization%option, &
                                                      surf_region, &
                                                      surf_region%filename)
#endif      
        endif
      else if (index(surf_region%filename,'.ss') > 0) then
        surf_region%sideset => RegionCreateSideset()
        call RegionReadFromFile(surf_region%sideset,surf_region%filename, &
                                surf_realization%option)
      else
        call RegionReadFromFile(surf_region,surf_realization%option, &
                                surf_region%filename)
      endif
    endif
    surf_region => surf_region%next
  enddo

end subroutine SurfaceInitReadRegionFiles

end module Surface_Init_module

#endif
