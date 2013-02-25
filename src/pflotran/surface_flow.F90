#ifdef SURFACE_FLOW

module Surface_Flow_module

  use Global_Aux_module
  !use Surface_Flow_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceFlowSetup, &
         SurfaceFlowTimeCut, &
         SurfaceFlowInitializeTimestep, &
         SurfaceFlowReadRequiredCardsFromInput, &
         SurfaceFlowRead, &
         SurfaceFlowKinematic, &
         SurfaceFlowKinematicDerivative, &
         SurfaceFlowDiffusion, &
         SurfaceFlowDiffusionDerivative, &
         SurfaceFlowResidual, &
         SurfaceFlowJacobian, &
         SurfaceFlowMaxChange, &
         SurfaceFlowUpdateSolution, &
         SurfaceFlowRHSFunction, &
         SurfaceFlowComputeMaxDt, &
         SurfaceFlowGetTecplotHeader
  
contains

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowSetup(surf_realization)

  use Surface_Realization_class
  
  type(surface_realization_type) :: surf_realization

  call SurfaceFlowSetPlotVariables(surf_realization)
  
end subroutine SurfaceFlowSetup

! ************************************************************************** !
!> This routine adds variables to be printed to list
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/30/12
! ************************************************************************** !
subroutine SurfaceFlowSetPlotVariables(surf_realization)
  
  use Surface_Realization_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => surf_realization%output_option%output_variable_list
  
  name = 'P'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_FLOW_PRESSURE)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceFlowSetPlotVariables

! ************************************************************************** !
!> This routine reads required surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceFlowReadRequiredCardsFromInput(surf_realization,input,option)

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
  !if (InputCheckExit(input,option)) exit    
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

end subroutine SurfaceFlowReadRequiredCardsFromInput

! ************************************************************************** !
!> This routine reads surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceFlowRead(surf_realization,surf_flow_solver,input,option)

  use Option_module
  use Input_module
  use String_module
  use Surface_Material_module
  use Surface_Realization_class
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Dataset_Aux_module
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
  type(dataset_type), pointer                  :: dataset

  type(patch_type), pointer                    :: patch
  type(output_option_type), pointer            :: output_option
  PetscReal :: units_conversion

  character(len=MAXWORDLENGTH)                 :: word
  character(len=MAXWORDLENGTH)                 :: card
  character(len=1) :: backslash

  PetscBool :: velocities
  PetscBool :: fluxes
  PetscBool :: continuation_flag
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
      dataset => DatasetCreate()
      call InputReadWord(input,option,dataset%name,PETSC_TRUE)
      call InputDefaultMsg(input,option,'Dataset name')
      call DatasetRead(dataset,input,option)
      call DatasetAddToList(dataset,surf_realization%datasets)
      nullify(dataset)
      
      !.........................................................................
      case default
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select
  enddo

  !surf_flow_solver%print_detailed_convergence = PETSC_TRUE

end subroutine SurfaceFlowRead


! ************************************************************************** !
!> This routine computes the residual equation.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowResidual(snes,xx,r,surf_realization,ierr)

  use Surface_Realization_class
  use Surface_Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(surface_realization_type) :: surf_realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string,string2
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
  
  surf_field => surf_realization%surf_field
  discretization => surf_realization%discretization
  option => surf_realization%option

  call DiscretizationGlobalToLocal(discretization,xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowResidualPatch1(snes,xx,r,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowResidualPatch2(snes,xx,r,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  surf_realization%iter_count = surf_realization%iter_count+1
  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

  if (surf_realization%debug%vecview_residual) then
    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rresidual_' // trim(adjustl(string2)) // '.out'
    !call PetscViewerASCIIOpen(surf_realization%option%mycomm,string, &
    !                          viewer,ierr)
    !call VecView(r,viewer,ierr)
    !call PetscViewerDestroy(viewer,ierr)

    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rresidual_' // trim(adjustl(string2)) // '.bin'
    write(*,*),'Writing -- ',trim(adjustl(string))
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    !call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)
  endif

  if (surf_realization%debug%vecview_solution) then
    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rxx_' // trim(adjustl(string2)) // '.out'
    !call PetscViewerASCIIOpen(surf_realization%option%mycomm,string, &
    !                          viewer,ierr)
    !call VecView(xx,viewer,ierr)
    !call PetscViewerDestroy(viewer,ierr)

    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rxx_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine SurfaceFlowResidual

! ************************************************************************** !
!> This routine computes the interior flux and boundary flux terms of the 
!! residual equation.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowResidualPatch1(snes,xx,r,surf_realization,ierr)

  use Water_EOS_module

  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(surface_realization_type) :: surf_realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), mannings_loc_p(:),xx_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: dP


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: slope, slope_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscViewer :: viewer
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  PetscReal :: dx, dy, dz, dist
  PetscReal :: vel
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

  r_p = 0.d0  

  call density(option%reference_temperature,option%reference_pressure,rho)
  !call nacl_den(option%reference_temperature,option%reference_pressure*1d-6,0.d0,rho)
  !rho = rho * 1.d3
  
  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      hw_up = xx_loc_p(ghosted_id_up)
      hw_dn = xx_loc_p(ghosted_id_dn)
      
      if(hw_up<0) then
        hw_up = 0.d0
        xx_loc_p(ghosted_id_up) = 0.d0
      endif
      if(hw_dn<0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
      endif
      
      if (hw_up<0 .or. hw_dn<0) then
        option%io_buffer = 'Surface water head negative'
        call printErrMsg(option)        
      endif
      
      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          call SurfaceFlowKinematic(hw_up,mannings_loc_p(ghosted_id_up), &
                                    hw_dn,mannings_loc_p(ghosted_id_dn), &
                                    slope, cur_connection_set%area(iconn), &
                                    option,vel,Res)
        case (DIFFUSION_WAVE)
          call SurfaceFlowDiffusion(hw_up,zc(ghosted_id_up), &
                                    mannings_loc_p(ghosted_id_up), &
                                    hw_dn,zc(ghosted_id_dn), &
                                    mannings_loc_p(ghosted_id_dn), &
                                    dist, &
                                    cur_connection_set%area(iconn), &
                                    option,vel,Res)
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(sum_connection) = Res(1)
      
      if (local_id_up>0) then
        r_p(local_id_up) = r_p(local_id_up) + Res(1)
      endif
         
      if (local_id_dn>0) then
        r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      hw_dn = xx_loc_p(ghosted_id_dn)
      if(hw_dn<0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
        write(*,*),'setting pressure values to zero for >> ',ghosted_id_dn
      endif

      call SurfaceBCFlux( boundary_condition%flow_condition%itype, &
                          hw_dn,slope_dn,mannings_loc_p(ghosted_id_dn), &
                          cur_connection_set%area(iconn),option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(sum_connection) = Res(1)
      
      r_p(local_id) = r_p(local_id) - Res(1)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

end subroutine SurfaceFlowResidualPatch1  

! ************************************************************************** !
!> This routine computes the accumulation and source/sink terms of the 
!! residual equation.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowResidualPatch2(snes,xx,r,surf_realization,ierr)

  use Water_EOS_module

  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use String_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(surface_realization_type) :: surf_realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), mannings_loc_p(:)
  PetscReal, pointer :: accum_p(:),area_p(:),xx_loc_p(:),qsrc_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: slope_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscViewer :: viewer
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: head         ! water height [m]
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  PetscReal :: dx, dy, dz
  PetscReal :: qsrc, qsrc_flow
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  
  call GridVecGetArrayF90(grid,r,r_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%flow_accum,accum_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)
  
  call density(option%reference_temperature,option%reference_pressure,rho)

  ! Accumulation terms ------------------------------------
  r_p = r_p - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle


    head = xx_loc_p(ghosted_id)
    if(head<0) head = 0.d0
    call SurfaceFlowAccumulation(head,area_p(local_id),option,Res)
    
    r_p(local_id) = r_p(local_id) + Res(1)
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    qsrc_flow = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc = m^3/sec
          qsrc = qsrc_flow*area_p(local_id)
        case(DISTRIBUTED_VOLUMETRIC_RATE_SS)
          ! qsrc = m^3/sec
          qsrc = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)*area_p(local_id)
        case default
          option%io_buffer = 'Source/Sink flow condition type not recognized'
          call printErrMsg(option)
      end select
      
      r_p(local_id) = r_p(local_id) - qsrc

    enddo
    source_sink => source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,r,r_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%flow_accum,accum_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)

  end subroutine SurfaceFlowResidualPatch2

! ************************************************************************** !
!> This routine computes the Jacbian.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowJacobian(snes,xx,A,B,flag,surf_realization,ierr)

  use Surface_Realization_class
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(surface_realization_type) :: surf_realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  character(len=MAXSTRINGLENGTH) :: string,string2

  call PetscLogEventBegin(logging%event_r_jacobian,ierr)

  option => surf_realization%option

  flag = SAME_NONZERO_PATTERN
  call MatGetType(A,mat_type,ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr)

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch

      call SurfaceFlowJacobianPatch1(snes,xx,J,J,flag,surf_realization,ierr)
      call SurfaceFlowJacobianPatch2(snes,xx,J,J,flag,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

  if (surf_realization%debug%matview_Jacobian) then
    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rjacobian_' // trim(adjustl(string2)) // '.out'
    !call PetscViewerASCIIOpen(surf_realization%option%mycomm,string, &
    !                          viewer,ierr)
    !call MatView(J,viewer,ierr)
    !call PetscViewerDestroy(viewer,ierr)

    !write(string,*) surf_realization%iter_count
    string = 'Surf_Rjacobian_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine SurfaceFlowJacobian


! ************************************************************************** !
!> This routine computes the interior flux and boundary flux terms of the
!! jacobian.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !

subroutine SurfaceFlowJacobianPatch1(snes,xx,A,B,flag,surf_realization,ierr)
       
  use Water_EOS_module

  use Connection_module
  use Surface_Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Surface_Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(surface_realization_type) :: surf_realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscInt :: icap_up,icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(surf_realization%option%nflowdof,surf_realization%option%nflowdof), &
               Jdn(surf_realization%option%nflowdof,surf_realization%option%nflowdof)
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  PetscReal, pointer :: mannings_loc_p(:),xx_loc_p(:)
  PetscReal :: dx, dy, dz, dist
  PetscReal :: slope_dn, slope
  
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: dP
  
  PetscViewer :: viewer

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option

  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

  call density(option%reference_temperature,option%reference_pressure,rho)

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------  
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
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist


      hw_up = xx_loc_p(ghosted_id_up)
      hw_dn = xx_loc_p(ghosted_id_dn)
      if(hw_up<0) hw_up = 0.d0
      if(hw_dn<0) hw_dn = 0.d0

      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          call SurfaceFlowKinematicDerivative(hw_up, mannings_loc_p(ghosted_id_up), &
                                              hw_dn, mannings_loc_p(ghosted_id_dn), &
                                              slope, &
                                              cur_connection_set%area(iconn), &
                                              option,Jup,Jdn)
        case (DIFFUSION_WAVE)
          call SurfaceFlowDiffusionDerivative(hw_up, &
                                              zc(ghosted_id_up), &
                                              mannings_loc_p(ghosted_id_up), &
                                              hw_dn, &
                                              zc(ghosted_id_dn), &
                                              mannings_loc_p(ghosted_id_dn), &
                                              dist, &
                                              cur_connection_set%area(iconn), &
                                              option,Jup,Jdn)
      end select
      if (local_id_up > 0) then
        call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                               Jup,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               Jdn,ADD_VALUES,ierr)
      endif

      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                              Jdn,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               Jup,ADD_VALUES,ierr)
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next

  enddo
   
  ! Boundary Flux Terms -----------------------------------  
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
  
      dx = xc(ghosted_id) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)


      hw_dn = xx_loc_p(ghosted_id)
      call SurfaceBCFluxDerivative( boundary_condition%flow_condition%itype, &
                                    hw_dn,slope_dn,mannings_loc_p(ghosted_id), &
                                    cur_connection_set%area(iconn),option,Jdn)

      call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                             ADD_VALUES,ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo


  call GridVecRestoreArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

end subroutine SurfaceFlowJacobianPatch1

! ************************************************************************** !
!> This routine computes the accumulation and source/sink terms of the 
!! jacobian.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !

subroutine SurfaceFlowJacobianPatch2(snes,xx,A,B,flag,surf_realization,ierr)
       
  use Water_EOS_module

  use Connection_module
  use Surface_Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Surface_Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(surface_realization_type) :: surf_realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscInt :: icap_up,icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(surf_realization%option%nflowdof,surf_realization%option%nflowdof), &
               Jdn(surf_realization%option%nflowdof,surf_realization%option%nflowdof)
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: head, hw_dn ! water height [m]
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  PetscReal, pointer :: mannings_loc_p(:),xx_loc_p(:),area_p(:)
  PetscReal :: dx, dy, dz
  PetscReal :: slope_dn, slope
  
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  
  PetscViewer :: viewer

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option

  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

  call density(option%reference_temperature,option%reference_pressure,rho)

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle


    head = xx_loc_p(ghosted_id)
    call SurfaceFlowAccumulationDerivative(head,area_p(local_id),option,Jup(1,1))
    
    call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                           ADD_VALUES,ierr)
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0

  call GridVecRestoreArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

end subroutine SurfaceFlowJacobianPatch2

! ************************************************************************** !
!> This routine computes the internal flux term for the residual under
!! kinematic-wave assumption.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowKinematic(hw_up, &
                                mannings_up, &
                                hw_dn, &
                                mannings_dn, &
                                slope, &
                                length, &
                                option, &
                                vel, &
                                Res)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: slope
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s
  
  PetscReal :: flux       ! units: m^2/s
  
  ! initialize
  flux = 0.d0
  vel  = 0.d0
  
  if (slope<0.d0) then
    vel =  sqrt(dabs(slope))/mannings_up*((hw_up)**(2.d0/3.d0))
    flux=  hw_up*vel
  else
    vel = -sqrt(dabs(slope))/mannings_dn*((hw_dn)**(2.d0/3.d0))
    flux=  hw_dn*vel
  endif

  Res(1) = flux*length

end subroutine SurfaceFlowKinematic

! ************************************************************************** !
!> This routine computes the derivative of the internal flux term for the
!! Jacobian under kinematic-wave assumption.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowKinematicDerivative(hw_up,mannings_up, &
                                          hw_dn,mannings_dn, &
                                          slope,length,option,Jup,Jdn)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: slope_up, slope_dn,slope
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: length
  PetscReal :: Jup(option%nflowdof,option%nflowdof), &
               Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: flux_dh_up, flux_dh_dn

  flux_dh_up = 0.d0
  flux_dh_dn = 0.d0
  
  if (slope < 0.d0) then
    flux_dh_up =  5.d0/3.d0 * sqrt(dabs(slope))/mannings_up*((hw_up)**(2.d0/3.d0))
  else
    flux_dh_dn = -5.d0/3.d0 * sqrt(dabs(slope))/mannings_dn*((hw_dn)**(2.d0/3.d0))
  endif
  
  Jup = flux_dh_up*length
  Jdn = flux_dh_dn*length

end subroutine SurfaceFlowKinematicDerivative

! ************************************************************************** !
!> This routine computes the internal flux term for the residual under
!! diffusion-wave assumption.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceFlowDiffusion(hw_up, &
                                zc_up, &
                                mannings_up, &
                                hw_dn, &
                                zc_dn, &
                                mannings_dn, &
                                dist, &
                                length, &
                                option, &
                                vel, &
                                Res)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: head_up, head_dn
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s

  PetscReal :: flux       ! units: m^2/s
  PetscReal :: Cd
  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: dhead

  ! initialize
  flux = 0.d0
  Cd = 1.0d0

  head_up = hw_up + zc_up
  head_dn = hw_dn + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (hw_up>0.d0) then
      !if (hw_up*Cd>hw_dn) then
      !  hw_half = 0.5d0*(hw_up+hw_dn)
      !else
        hw_half = hw_up
      !endif
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    if (hw_dn>0.d0) then
      !if (hw_dn*Cd>hw_up) then
      !  hw_half = 0.5d0*(hw_up+hw_dn)
      !else
        hw_half = hw_dn
      !endif
    else
      hw_half = 0.d0
    endif
  endif
  
  !vel = -dsign(1.d0,head_dn-head_up)*(1.0d0/mannings_half)* &
  !          (hw_half**(2.d0/3.d0))* &
  !          (abs((head_dn-head_up)/dist)**(1.d0/2.d0))
  dhead=head_up-head_dn
  if(abs(dhead)<eps) then
    dhead=0.d0
    vel = 0.d0
  else
!    vel = (dhead)/mannings_half/(dist**(1.d0/2.d0))* &
!          (hw_half**(2.d0/3.d0))* &
!          1.d0/(abs(dhead)**(1.d0/2.d0))
    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)
  endif

  flux = hw_half*vel
  Res(1) = flux*length

end subroutine SurfaceFlowDiffusion

! ************************************************************************** !
!> This routine computes the derivative of the internal flux term for the
!! Jacobian under diffusion-wave assumption.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/06/12
! ************************************************************************** !
subroutine SurfaceFlowDiffusionDerivative(hw_up,zc_up,mannings_up, &
                                          hw_dn,zc_dn,mannings_dn, &
                                          dist,length,option,Jup,Jdn)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: dist,length
  PetscReal :: Jup(option%nflowdof,option%nflowdof), &
               Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: flux_dh_up, flux_dh_dn
  PetscReal :: Cd
  PetscReal :: mannings_half,hw_half
  PetscReal :: dhw_half_dhw_dn,dhw_half_dhw_up
  PetscReal :: term0,term1,term2,term3,term4
  PetscReal :: head_up, head_dn, dhead
  PetscReal :: Jup_new(option%nflowdof,option%nflowdof), &
               Jdn_new(option%nflowdof,option%nflowdof)
  PetscReal :: res_pert_up(1),res_pert_dn(1),vel,Jup_pert,Jdn_pert
  PetscReal :: res(1:option%nflowdof)   ! units: m^3/s
  PetscReal :: dhead_dn_dhw_dn, dhead_up_dhw_up

  flux_dh_up = 0.d0
  flux_dh_dn = 0.d0
  
  ! initialize
  Cd = 1.0d0
  dhw_half_dhw_up = 0.d0
  dhw_half_dhw_dn = 0.d0

  head_up = hw_up + zc_up
  head_dn = hw_dn + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (hw_up>0.d0) then
      !if (hw_up*Cd>hw_dn) then
      !  hw_half = 0.5d0*(hw_up+hw_dn)
      !  dhw_half_dhw_dn = 0.5d0
      !  dhw_half_dhw_up = 0.5d0
      !else
        hw_half         = hw_up
        dhw_half_dhw_up = 1.d0
      !endif
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    if (hw_dn>0.d0) then
      !if (hw_dn*Cd>hw_up) then
      !  hw_half = 0.5d0*(hw_up+hw_dn)
      !  dhw_half_dhw_up = 0.5d0
      !  dhw_half_dhw_dn = 0.5d0
      !else
        hw_half = hw_dn
        dhw_half_dhw_dn = 1.d0
      !endif
    else
      hw_half = 0.d0
    endif
  endif

  dhead = head_up-head_dn

  if(abs(dhead)<eps) then
    dhw_half_dhw_up=0.d0
    dhw_half_dhw_dn=0.d0
    dhead=0.d0
    term0=0.d0
    Jup=0.d0
    Jdn=0.d0
    return
  endif


  term0 = (2.0d0*dhead**2.d0 - head_dn**2.d0 - head_up**2.d0 + 2*head_dn*head_up)
  term0 = term0/(2.d0*((dhead)**2.d0)*(abs(dhead)**(1.d0/2.d0)))

  term1 = (5.d0/3.d0)*(hw_half**(2.d0/3.d0))/mannings_half/(dist**(1.d0/2.d0))
  term1 = term1*dhead/(abs(dhead)**(1.d0/2.d0))*dhw_half_dhw_up

  term2 = (hw_half**(5.d0/3.d0))/mannings_half/(dist**(1.d0/2.d0))
  term2 = term0*term2

  Jup = (term1+term2)*length

  term1 = (5.d0/3.d0)*(hw_half**(2.d0/3.d0))/mannings_half/(dist**(1.d0/2.d0))
  term1 = term1*dhead/(abs(dhead)**(1.d0/2.d0))*dhw_half_dhw_dn

  term2 = (hw_half**(5.d0/3.d0))/mannings_half/(dist**(1.d0/2.d0))
  term2 = -term0*term2

  Jdn = (term1+term2)*length

  if(option%numerical_derivatives_flow) then
    call SurfaceFlowDiffusion(hw_up,zc_up,mannings_up, &
                              hw_dn,zc_dn,mannings_dn, &
                              dist,length,option,vel,res)

    call SurfaceFlowDiffusion(hw_up+perturbation_tolerance,zc_up,mannings_up, &
                              hw_dn,zc_dn,mannings_dn, &
                              dist,length,option,vel,res_pert_up)

    call SurfaceFlowDiffusion(hw_up,zc_up,mannings_up, &
                              hw_dn+perturbation_tolerance,zc_dn,mannings_dn, &
                              dist,length,option,vel,res_pert_dn)
    Jup_pert=(res_pert_up(1)-res(1))/perturbation_tolerance
    Jdn_pert=(res_pert_dn(1)-res(1))/perturbation_tolerance

    Jup=Jup_pert
    Jdn=Jdn_pert
  endif


end subroutine SurfaceFlowDiffusionDerivative

! ************************************************************************** !
!> This routine resets arrays for time step cut.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/31/12
! ************************************************************************** !
subroutine SurfaceFlowTimeCut(surf_realization)
 
  use Surface_Realization_class
  use Surface_Field_module
 
  implicit none
  
  type(surface_realization_type) :: surf_realization
  type(surface_field_type), pointer :: surf_field
  
  PetscErrorCode :: ierr

  surf_field => surf_realization%surf_field

  call VecCopy(surf_field%flow_yy,surf_field%flow_xx,ierr)
  call SurfaceFlowInitializeTimestep(surf_realization)
 
end subroutine SurfaceFlowTimeCut

! ************************************************************************** !
!> This routine updates the data prior to time step.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowInitializeTimestep(surf_realization)

  use Surface_Realization_class
  use Surface_Field_module
  
  implicit none

  type(surface_realization_type) :: surf_realization

  call SurfaceFlowUpdateFixedAccum(surf_realization)

end subroutine SurfaceFlowInitializeTimestep

! ************************************************************************** !
!> This routine updates the fixed portion of the accumulation term
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateFixedAccum(surf_realization)

  use Surface_Realization_class
  use Level_module
  use Patch_module

  type(surface_realization_type) :: surf_realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowUpdateFixedAccumPatch(surf_realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceFlowUpdateFixedAccum

! ************************************************************************** !
!> This routine updates the fixed portion of the accumulation term over a 
!! single patch
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateFixedAccumPatch(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Grid_module
  use Surface_Field_module
  use Water_EOS_module

  implicit none

  type(surface_realization_type) :: surf_realization

  type(option_type), pointer         :: option
  type(patch_type), pointer          :: patch
  type(grid_type), pointer           :: grid
  type(surface_field_type),pointer   :: surf_field
  
  PetscInt             :: local_id, ghosted_id
  PetscReal, pointer   :: accum_p(:),area_p(:),xx_p(:)
  PetscReal            :: rho          ! density      [kg/m^3]
  PetscReal            :: head         ! [m]

  PetscErrorCode :: ierr

  option     => surf_realization%option
  grid       => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid, surf_field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%area, area_p,ierr)
  call GridVecGetArrayF90(grid, surf_field%flow_xx, xx_p, ierr)

  call density(option%reference_temperature,option%reference_pressure,rho)

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)
    
    head = xx_p(local_id)
    call SurfaceFlowAccumulation(head,area_p(local_id),option, &
                                 accum_p(local_id:local_id))
  enddo

  call GridVecRestoreArrayF90(grid ,surf_field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%area, area_p,ierr)
  call GridVecRestoreArrayF90(grid, surf_field%flow_xx, xx_p, ierr)

end subroutine SurfaceFlowUpdateFixedAccumPatch

! ************************************************************************** !
!> This routine computes the non-fixed portion of the accumulation term
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowAccumulation(head, area, option, Res)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal         :: Res(1:option%nflowdof)
  PetscReal         :: head, area

  Res(1) = head*area/option%surf_flow_dt

end subroutine SurfaceFlowAccumulation

! ************************************************************************** !
!> This routine computes the derivative of the non-fixed portion of the 
!! accumulation term for the Jacobian
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowAccumulationDerivative(head, area, option, J)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal         :: J(option%nflowdof,option%nflowdof)
  PetscReal         :: head, area

  J(1,1) = area /option%surf_flow_dt

end subroutine SurfaceFlowAccumulationDerivative


! ************************************************************************** !
!> This routine computes the boundary flux term for the residual
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceBCFlux(ibndtype,head,slope,mannings, &
                         length,option,vel,Res)

  use Option_module
  
  implicit none

  type(option_type) :: option
  PetscReal :: head
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt  :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: Res(1:option%nflowdof) 

  PetscInt :: pressure_bc_type

  flux = 0.d0
  vel = 0.d0
  
  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
      else
        vel = -sqrt(dabs(slope))/mannings*((head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = head*vel
  Res(1) = flux*length

end subroutine SurfaceBCFlux

! ************************************************************************** !
!> This routine computes the derivative of the boundary flux term for the 
!! Jacobian
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceBCFluxDerivative(ibndtype,head,slope,mannings, &
                                   length,option,J)

  use Option_module
  
  implicit none

!#include "definitions.h"

  type(option_type) :: option
  PetscReal :: J(1:option%nflowdof)
  PetscReal :: head
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscInt :: ibndtype(:)

  PetscInt :: pressure_bc_type
  PetscReal :: flux_dh
  PetscReal :: J_pert(1:option%nflowdof),vel
  PetscReal :: res(1)
  PetscReal :: res_pert(1)

  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        flux_dh = 0.d0
      else
        flux_dh = -5.d0/3.d0*sqrt(dabs(slope))/mannings*((head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select

  J(1) = flux_dh*length

  if(option%numerical_derivatives_flow) then
    call SurfaceBCFlux(ibndtype,head,slope,mannings, &
                       length,option,vel,res)

    call SurfaceBCFlux(ibndtype,head+perturbation_tolerance,slope,mannings, &
                       length,option,vel,res_pert)

    J_pert(1)=(res_pert(1)-res(1))/perturbation_tolerance
    J(1)=J_pert(1)
  endif


end subroutine SurfaceBCFluxDerivative

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceFlowRHSFunction(ts,t,xx,ff,surf_realization,ierr)

  use Surface_Realization_class
  use Surface_Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  PetscViewer :: viewer

  type(discretization_type), pointer   :: discretization
  type(surface_field_type), pointer    :: surf_field
  type(level_type), pointer            :: cur_level
  type(patch_type), pointer            :: cur_patch
  type(option_type), pointer           :: option
  character(len=MAXSTRINGLENGTH)       :: string,string2

  write(*,*),'In SurfaceFlowRHSFunction:           ',t

  surf_realization%iter_count = surf_realization%iter_count+1
  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

!  call VecView(ff,PETSC_VIEWER_STDOUT_WORLD,ierr)

  surf_field      => surf_realization%surf_field
  discretization  => surf_realization%discretization
  option          => surf_realization%option

  call DiscretizationGlobalToLocal(discretization,xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowRHSFunctionPatch1(ts,t,xx,ff,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowRHSFunctionPatch2(ts,t,xx,ff,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

!  call VecView(ff,PETSC_VIEWER_STDOUT_WORLD,ierr)
  if (surf_realization%debug%vecview_solution) then
    string = 'Surf_xx_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    string = 'Surf_ff_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(ff,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  endif

end subroutine SurfaceFlowRHSFunction

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceFlowRHSFunctionPatch1(ts,t,xx,ff,surf_realization,ierr)

  use water_eos_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: max_allowable_dt

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),xx_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid,ff,ff_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%mannings_loc,mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

  ff_p = 0.d0
  Res  = 0.d0
  max_allowable_dt = 1.d10

  call density(option%reference_temperature,option%reference_pressure,rho)
  !call nacl_den(option%reference_temperature,option%reference_pressure*1d-6,0.d0,rho)
  !rho = rho * 1.d3

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      hw_up = xx_loc_p(ghosted_id_up)
      hw_dn = xx_loc_p(ghosted_id_dn)
      
      if(hw_up<0.d0) then
        hw_up = 0.d0
        xx_loc_p(ghosted_id_up) = 0.d0
      endif
      if(hw_dn<0.d0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
      endif
      
      if (hw_up<0.d0 .or. hw_dn<0.d0) then
        option%io_buffer = 'Surface water head negative'
        call printErrMsg(option)        
      endif

      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          !call SurfaceFlowKinematic(hw_up,mannings_loc_p(ghosted_id_up), &
          !                          hw_dn,mannings_loc_p(ghosted_id_dn), &
          !                          slope, cur_connection_set%area(iconn), &
          !                          option,vel,Res)
        case (DIFFUSION_WAVE)
          call SurfaceFlowDiffusion(hw_up,zc(ghosted_id_up), &
                                    mannings_loc_p(ghosted_id_up), &
                                    hw_dn,zc(ghosted_id_dn), &
                                    mannings_loc_p(ghosted_id_dn), &
                                    dist, &
                                    cur_connection_set%area(iconn), &
                                    option,vel,Res)
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(sum_connection) = Res(1)
      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)
      
      if (local_id_up>0) then
        ff_p(local_id_up) = ff_p(local_id_up) - Res(1)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        ff_p(local_id_dn) = ff_p(local_id_dn) + Res(1)/area_p(local_id_dn)
      endif


    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      hw_dn = xx_loc_p(ghosted_id_dn)
      if(hw_dn<0.d0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
        write(*,*),'setting pressure values to zero for >> ',ghosted_id_dn
      endif

      call SurfaceBCFlux( boundary_condition%flow_condition%itype, &
                          hw_dn,slope_dn,mannings_loc_p(ghosted_id_dn), &
                          cur_connection_set%area(iconn),option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(sum_connection) = Res(1)
      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)
      
      ff_p(local_id) = ff_p(local_id) + Res(1)/area_p(local_id)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,ff,ff_p, ierr)
  call GridVecRestoreArrayF90(grid,surf_field%mannings_loc,mannings_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)

end subroutine SurfaceFlowRHSFunctionPatch1

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceFlowRHSFunctionPatch2(ts,t,xx,ff,surf_realization,ierr)

  use water_eos_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                 :: grid
  type(patch_type), pointer                :: patch
  type(option_type), pointer               :: option
  type(surface_field_type), pointer        :: surf_field
  type(coupler_type), pointer              :: source_sink
  type(connection_set_list_type), pointer  :: connection_set_list
  type(connection_set_type), pointer       :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt :: local_id, ghosted_id

  PetscReal :: qsrc, qsrc_flow
  PetscReal, pointer ::ff_p(:),area_p(:),xx_loc_p(:)

  patch      => surf_realization%patch
  grid       => patch%grid
  option     => surf_realization%option
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid,ff,ff_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    if(source_sink%flow_condition%rate%itype/=DISTRIBUTED_VOLUMETRIC_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=DISTRIBUTED_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc = m^3/sec
          qsrc = qsrc_flow*area_p(local_id)
        case(DISTRIBUTED_VOLUMETRIC_RATE_SS)
          ! qsrc = m^3/sec
          qsrc = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)*area_p(local_id)
        case default
          option%io_buffer = 'Source/Sink flow condition type not recognized'
          call printErrMsg(option)
      end select
      
      ff_p(local_id) = ff_p(local_id) + qsrc/area_p(local_id)

    enddo
    source_sink => source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,ff,ff_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)

end subroutine SurfaceFlowRHSFunctionPatch2

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceFlowComputeMaxDt(surf_realization,max_allowable_dt)

  use water_eos_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: max_allowable_dt

  PetscReal, pointer :: mannings_loc_p(:),xx_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid,surf_field%mannings_loc,mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10

  call density(option%reference_temperature,option%reference_pressure,rho)
  !call nacl_den(option%reference_temperature,option%reference_pressure*1d-6,0.d0,rho)
  !rho = rho * 1.d3

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      hw_up = xx_loc_p(ghosted_id_up)
      hw_dn = xx_loc_p(ghosted_id_dn)
      
      if(hw_up<0.d0) then
        hw_up = 0.d0
        xx_loc_p(ghosted_id_up) = 0.d0
      endif
      if(hw_dn<0.d0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
      endif
      
      if (hw_up<0.d0 .or. hw_dn<0.d0) then
        option%io_buffer = 'Surface water head negative'
        call printErrMsg(option)        
      endif

      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          !call SurfaceFlowKinematic(hw_up,mannings_loc_p(ghosted_id_up), &
          !                          hw_dn,mannings_loc_p(ghosted_id_dn), &
          !                          slope, cur_connection_set%area(iconn), &
          !                          option,vel,Res)
        case (DIFFUSION_WAVE)
          call SurfaceFlowDiffusion(hw_up,zc(ghosted_id_up), &
                                    mannings_loc_p(ghosted_id_up), &
                                    hw_dn,zc(ghosted_id_dn), &
                                    mannings_loc_p(ghosted_id_dn), &
                                    dist, &
                                    cur_connection_set%area(iconn), &
                                    option,vel,Res)
      end select

      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      hw_dn = xx_loc_p(ghosted_id_dn)
      if(hw_dn<0.d0) then
        hw_dn = 0.d0
        xx_loc_p(ghosted_id_dn) = 0.d0
        write(*,*),'setting pressure values to zero for >> ',ghosted_id_dn
      endif

      call SurfaceBCFlux( boundary_condition%flow_condition%itype, &
                          hw_dn,slope_dn,mannings_loc_p(ghosted_id_dn), &
                          cur_connection_set%area(iconn),option,vel,Res)

      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,surf_field%mannings_loc,mannings_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)

end subroutine SurfaceFlowComputeMaxDt

! ************************************************************************** !
!> This routine computes the maximum change in the solution vector
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/21/12
! ************************************************************************** !
subroutine SurfaceFlowMaxChange(surf_realization)

  use Surface_Realization_class
  use Option_module
  use Surface_Field_module
  
  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(surface_field_type), pointer  :: surf_field  
  
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  call VecWAXPY(surf_field%flow_dxx,-1.d0,surf_field%flow_xx,surf_field%flow_yy,ierr)
  call VecStrideNorm(surf_field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

end subroutine SurfaceFlowMaxChange

! ************************************************************************** !
!> This routine updates data in module after a successful time step
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/22/12
! ************************************************************************** !
subroutine SurfaceFlowUpdateSolution(surf_realization)

  use Surface_Realization_class
  use Surface_Field_module

  implicit none

  type(surface_realization_type)   :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode                   :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr)

end subroutine SurfaceFlowUpdateSolution

! ************************************************************************** !
!> This routine surface flow tecplot file header
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/29/12
! ************************************************************************** !
function SurfaceFlowGetTecplotHeader(surf_realization,icolumn)

  use Surface_Realization_class
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: SurfaceFlowGetTecplotHeader
  type(surface_realization_type) :: surf_realization
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: string, string2

  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [m]"'')')
  endif
  string = trim(string) // trim(string2)
#ifdef GLENN_NEW_IO
  !call OutputOptionAddPlotVariable(realization%output_option,PRESSURE, &
  !                           ZERO_INTEGER,ZERO_INTEGER)
#endif

  SurfaceFlowGetTecplotHeader = string

end function SurfaceFlowGetTecplotHeader

end module Surface_Flow_module

#endif
