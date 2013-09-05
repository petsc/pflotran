#ifdef GEOMECH

module Geomechanics_Init_module

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"


  public :: GeomechicsInitReadRequiredCards, &
            GeomechanicsInitReadInput

contains

! ************************************************************************** !
!
! GeomechicsInitReadRequiredCards: Reads the required input file cards
! related to geomechanics
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechicsInitReadRequiredCards(geomech_realization)

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Option_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer          :: option
  type(input_type), pointer           :: input
  
  option         => geomech_realization%option  
  input          => geomech_realization%input
  
! Read in select required cards
!.........................................................................

  ! GEOMECHANICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if(InputError(input)) return
  option%ngeomechdof = 3  ! displacements in x, y, z directions
  
  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
  call GeomechanicsInit(geomech_realization,input,option)  


end subroutine GeomechicsInitReadRequiredCards

! ************************************************************************** !
!
! GeomechanicsInit: Reads the required geomechanics data from input file
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInit(geomech_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Grid_module
  
  implicit none
  
  type(geomech_realization_type)             :: geomech_realization
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_patch_type), pointer          :: patch
  type(input_type)                           :: input
  type(option_type), pointer                 :: option
  character(len=MAXWORDLENGTH)               :: word
  type(unstructured_grid_type), pointer      :: ugrid
  
  discretization       => geomech_realization%discretization
       
 input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  call InputReadFlotranString(input,option)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
  call StringToUpper(word)
    
  select case(trim(word))
    case ('TYPE')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword','TYPE')
      call StringToUpper(word)

      select case(trim(word))
        case ('UNSTRUCTURED')
          discretization%itype = UNSTRUCTURED_GRID
          call InputReadNChars(input,option, &
                               discretization%filename, &
                               MAXSTRINGLENGTH, &
                               PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','filename')

          discretization%grid  => GMGridCreate()
          ugrid => UGridCreate()
          call UGridRead(ugrid,discretization%filename,option)
          call UGridDecompose(ugrid,option)
          call CopySubsurfaceGridtoGeomechGrid(ugrid, &
                                               discretization%grid,option)
          patch => GeomechanicsPatchCreate()
          patch%geomech_grid => discretization%grid
          geomech_realization%geomech_patch => patch
        case default
          option%io_buffer = 'Geomechanics supports only unstructured grid'
          call printErrMsg(option)
      end select
  end select
     
end subroutine GeomechanicsInit

! ************************************************************************** !
!
! GeomechanicsInitReadInput: Reads the geomechanics input data 
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsInitReadInput(geomech_realization,geomech_solver, &
                                     input,option)

  use Option_module
  use Input_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Material_module
  use Geomechanics_Region_module
  use Geomechanics_Debug_module
  use Geomechanics_Strata_module
  use Geomechanics_Condition_module
  use Geomechanics_Coupler_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Solver_module
  use Units_module
  use Waypoint_module

  ! Still need to add other geomech modules for output, etc once created
  
  implicit none
  
  type(geomech_realization_type)               :: geomech_realization
  type(solver_type)                            :: geomech_solver
  type(input_type)                             :: input
  type(option_type)                            :: option
  
  type(geomech_discretization_type), pointer   :: discretization
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(waypoint_type), pointer                 :: waypoint
  type(geomech_grid_type), pointer             :: grid
  type(gm_region_type), pointer                :: region
  type(geomech_debug_type), pointer            :: debug
  type(geomech_strata_type), pointer           :: strata
  type(geomech_condition_type), pointer        :: condition
  type(geomech_coupler_type), pointer          :: coupler
  type(output_option_type), pointer            :: output_option
  PetscReal                                    :: units_conversion

  character(len=MAXWORDLENGTH)                 :: word
  character(len=MAXWORDLENGTH)                 :: card
  character(len=1)                             :: backslash
  
  PetscBool                                    :: continuation_flag
  PetscReal                                    :: temp_real, temp_real2
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''
  
  discretization => geomech_realization%discretization
  output_option => geomech_realization%output_option
    
  if (associated(geomech_realization%geomech_patch)) grid => &
    geomech_realization%geomech_patch%geomech_grid
    
  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    option%io_buffer = 'word :: ' // trim(word)
    call printMsg(option)   

    select case(trim(word))
    
      !.........................................................................
      ! Read geomechanics grid information
      case ('GEOMECHANICS_GRID')
        call InputSkipToEND(input,option,trim(word))
        
      !.........................................................................
      ! Read geomechanics material information
      case ('GEOMECHANICS_MATERIAL_PROPERTY')
        geomech_material_property => GeomechanicsMaterialPropertyCreate()

        call InputReadWord(input,option,geomech_material_property%name, &
                           PETSC_TRUE)
                           
        call InputErrorMsg(input,option,'name','GEOMECHANICS_MATERIAL_PROPERTY')
        call GeomechanicsMaterialPropertyRead(geomech_material_property,input, &
                                              option)
        call GeomechanicsMaterialPropertyAddToList(geomech_material_property, &
                                geomech_realization%geomech_material_properties)
        nullify(geomech_material_property)

      !.........................................................................
      case ('GEOMECHANICS_REGION')
        region => GeomechRegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_REGION')
        call printMsg(option,region%name)
        call GeomechRegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call GeomechRegionAddToList(region,geomech_realization%geomech_regions)
        nullify(region)
 
      !.........................................................................
      case ('GEOMECHANICS_CONDITION')
        condition => GeomechConditionCreate(option)
        call InputReadWord(input,option,condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'GEOMECHANICS_CONDITION','name')
        call printMsg(option,condition%name)
        call GeomechConditionRead(condition,input,option)
        call GeomechConditionAddToList(condition,geomech_realization%geomech_conditions)
        nullify(condition)
        
     !.........................................................................
      case ('GEOMECHANICS_BOUNDARY_CONDITION')
        coupler =>  GeomechCouplerCreate(GM_BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Geomech Boundary Condition name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
        
      !.........................................................................
      case ('GEOMECHANICS_SRC_SINK')
        coupler => GeomechCouplerCreate(GM_SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
                 
      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadNewton(geomech_solver,input,option)
        end select     
        
      !.........................................................................
      case ('GEOMECHANICS_DEBUG')
        call GeomechDebugRead(geomech_realization%debug,input,option)    

      !.........................................................................
      case('GEOMECHANICS_SUBSURFACE_COUPLING')
        option%geomech_subsurf_coupling = PETSC_TRUE        
        
      !.........................................................................
      case ('GEOMECHANICS_OUTPUT')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT')
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
              call InputErrorMsg(input,option,'units','GEOMECHANICS_OUTPUT')
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
                    write(*,*),'Inserting waypoint in geomech_realization: ',waypoint%time
                    call WaypointInsertInList(waypoint,geomech_realization%waypoints)
                  endif
                enddo
                if (.not.continuation_flag) exit
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
              enddo
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'GEOMECHANICS_OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC,OUTPUT_FILE')
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
                                     'GEOMECHANICS_OUTPUT,PERIODIC,SCREEN')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,SCREEN.'
                  call printErrMsg(option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'GEOMECHANICS_OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word,option)
                  output_option%periodic_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then

                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word,option)
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'GEOMECHANICS_OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_output = PETSC_TRUE
                        write(*,*),'Inserting waypoint in geomech_realization: >>>>>>>> ',waypoint%time
                        call WaypointInsertInList(waypoint,geomech_realization%waypoints)
                        temp_real = temp_real + output_option%periodic_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'GEOMECHANICS_OUTPUT,PERIODIC,TIME')
                    endif
                  endif
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in GEOMECHANICS_OUTPUT,PERIODIC,'// &
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
                                     'GEOMECHANICS_OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC_OBSERVATION,TIME')
                  units_conversion = UnitsConvertToInternal(word,option) 
                  output_option%periodic_tr_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_tr_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'GEOMECHANICS_OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,'// &
                                     'PERIODIC_OBSERVATION,TIMESTEP.'
                  call printErrMsg(option)
              end select
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT,FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'GEOMECHANICS_OUTPUT,FORMAT,HDF5,# FILES')
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
                  call InputErrorMsg(input,option,'TECPLOT','GEOMECHANICS_OUTPUT,FORMAT') 
                  call StringToUpper(word)
                  output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT ! By default it is unstructured
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,FORMAT.'
                  call printErrMsg(option)
              end select
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in GEOMECHANICS_OUTPUT.'
              call printErrMsg(option)
          end select
        enddo
                        
      !.........................................................................
      case ('GEOMECHANICS_STRATIGRAPHY','GEOMECHANICS_STRATA')
        strata => GeomechStrataCreate()
        call GeomechStrataRead(strata,input,option)
        call GeomechRealizAddStrata(geomech_realization,strata)
        nullify(strata)       
        
      !.........................................................................
      case default
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select
  enddo
  
end subroutine GeomechanicsInitReadInput

! ************************************************************************** !
!
! GeomechInitMatPropToGeomechRegions: This routine assigns geomech material 
!                                     properties to associated regions 
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechInitMatPropToGeomechRegions(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Discretization_module
  use Geomechanics_Strata_module
  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Option_module

  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: ivertex, local_id, ghosted_id, natural_id, geomech_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_field_type), pointer :: field
  type(geomech_strata_type), pointer :: strata
  type(geomech_patch_type), pointer :: patch  

  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: null_geomech_material_property
  type(gm_region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  
  option => geomech_realization%option
  discretization => geomech_realization%discretization
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch

  ! loop over all patches and allocation material id arrays
  if (.not.associated(patch%imat)) then
    allocate(patch%imat(patch%geomech_grid%ngmax_node))
    ! initialize to "unset"
    patch%imat = -999
  endif

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  grid => patch%geomech_grid
  strata => patch%geomech_strata%first
  do
    if (.not.associated(strata)) exit
    ! Read in cell by cell material ids if they exist
    if (.not.associated(strata%region) .and. strata%active) then
      option%io_buffer = 'Reading of material prop from file for' // &
        ' geomech is not implemented.'
      call printErrMsgByRank(option)
    ! Otherwise, set based on region
    else if (strata%active) then
      update_ghosted_material_ids = PETSC_TRUE
      region => strata%region
      geomech_material_property => strata%material_property
      if (associated(region)) then
        istart = 1
        iend = region%num_verts
      else
        istart = 1
        iend = grid%nlmax_node
      endif
      do ivertex = istart, iend
        if (associated(region)) then
          local_id = region%vertex_ids(ivertex)
        else
          local_id = ivertex
        endif
        ghosted_id = grid%nL2G(local_id)
        patch%imat(ghosted_id) = geomech_material_property%id
      enddo
    endif
    strata => strata%next
  enddo
    
  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call GeomechRealizLocalToLocalWithArray(geomech_realization, &
                                            MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_geomech_material_property => GeomechanicsMaterialPropertyCreate()
  do local_id = 1, grid%nlmax_node
    ghosted_id = grid%nL2G(local_id)
    geomech_material_id = patch%imat(ghosted_id)
    if (geomech_material_id == 0) then ! accomodate inactive cells
      geomech_material_property = null_geomech_material_property
    else if ( geomech_material_id > 0 .and. &
              geomech_material_id <= &
              size(geomech_realization%geomech_material_property_array)) then
      geomech_material_property => &
         geomech_realization% &
           geomech_material_property_array(geomech_material_id)%ptr
      if (.not.associated(geomech_material_property)) then
        write(dataset_name,*) geomech_material_id
        option%io_buffer = 'No material property for geomech material id ' // &
                            trim(adjustl(dataset_name)) &
                            //  ' defined in input file.'
        call printErrMsgByRank(option)
      endif
    else if (geomech_material_id < -998) then 
      write(dataset_name,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized geomech material id in patch at cell ' // &
                         trim(adjustl(dataset_name))
      call printErrMsgByRank(option)
    else if (geomech_material_id > size(geomech_realization% &
      geomech_material_property_array)) then
      write(option%io_buffer,*) geomech_material_id
      option%io_buffer = 'Unmatched geomech material id in patch:' // &
        adjustl(trim(option%io_buffer))
      call printErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with geomech material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn/Satish!'
      call printErrMsgByRank(option)
    endif
  enddo ! local_id - loop

  
  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)

end subroutine GeomechInitMatPropToGeomechRegions
 
end module Geomechanics_Init_module
#endif
! GEOMECH
