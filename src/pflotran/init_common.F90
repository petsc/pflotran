module Init_Common_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"

  public :: Init, &
            InitCommonReadRegionFiles, &
            InitCommonReadVelocityField, &
            InitCommonVerifyAllCouplers
            
contains

! ************************************************************************** !

subroutine Init(simulation)
  ! 
  ! Initializes pflotran
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Discretization_module
  use Realization_class
  use Material_module
  use Timestepper_module
  use Field_module
  use Connection_module   
  use Coupler_module
  use Debug_module
  use Convergence_module
  use Waypoint_module
  use Patch_module
! use Mass_Balance_module
  use Logging_module  
  use Reaction_Database_module
  use Reaction_Database_hpt_module
  use Input_Aux_module
  use Condition_Control_module
  
  use Flash2_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use Richards_module
  use TH_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  
  use Secondary_Continuum_module, only : SecondaryRTUpdateIterate
  
  use Global_module
  use Variables_module
  
  use EOS_module
  use EOS_Water_module
!  use Utility_module
  use Output_module
  use Output_Aux_module
  use Regression_module
    
  use Surface_Field_module
  use Surface_Flow_module
  use Surface_Global_module
  use Surface_Init_module !, only : SurfaceInitReadRequiredCards
  use Surface_Realization_class
  use Surface_TH_module
  use Grid_Unstructured_module

  use Geomechanics_Realization_class
  use Geomechanics_Init_module, only : GeomechicsInitReadRequiredCards
  use Geomechanics_Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Field_module
  use Geomechanics_Global_module
  use Geomechanics_Force_module

  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXSTRINGLENGTH) :: filename, filename_out

  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  type(realization_type), pointer :: realization
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(waypoint_list_type), pointer :: waypoint_list
  type(input_type), pointer :: input
  type(output_variable_type), pointer :: output_variable
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: temp_int
  PetscInt :: flowortranpc    
  PetscErrorCode :: ierr
  PCSide:: pcside
  PetscReal :: dum1
  PetscReal :: min_value
  SNESLineSearch :: linesearch
  type(timestepper_type), pointer               :: surf_flow_timestepper
  type(solver_type), pointer                :: surf_flow_solver
  type(surface_realization_type), pointer   :: surf_realization
  type(solver_type), pointer                :: geomech_solver
  type(timestepper_type), pointer               :: geomech_timestepper
  type(geomech_realization_type), pointer   :: geomech_realization

  ! set pointers to objects
  realization => simulation%realization
  discretization => realization%discretization
  option => realization%option
  input => realization%input
  surf_realization  => simulation%surf_realization
  surf_flow_timestepper => simulation%surf_flow_timestepper
  geomech_realization => simulation%geomech_realization
  geomech_timestepper => simulation%geomech_timestepper
  
  nullify(flow_solver)
  nullify(tran_solver)

  ! sets pointers to EOS procedures
  call EOSInit()
  
  if (OptionPrintToScreen(option)) then
    temp_int = 6
    call InitCommonPrintPFLOTRANHeader(option,temp_int)
  endif
  
  realization%input => InputCreate(IN_UNIT,option%input_filename,option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  if (OptionPrintToFile(option)) then
    call InitCommonPrintPFLOTRANHeader(option,option%fid_out)
  endif
  
  ! read required cards
  call InitSubsurfaceReadRequiredCards(realization)
  !geh: surf_realization%input is never freed
  surf_realization%input => InputCreate(IN_UNIT,option%input_filename,option)
  surf_realization%subsurf_filename = realization%discretization%filename
  call SurfaceInitReadRequiredCards(simulation%surf_realization)

  geomech_realization%input => InputCreate(IN_UNIT,option%input_filename,option)
  call GeomechicsInitReadRequiredCards(simulation%geomech_realization)

#if 0  
  patch => realization%patch

  if (associated(patch)) then
     if (associated(patch%grid)) then
        grid => patch%grid
     endif
  endif
#endif  

  ! process command line options
  call OptionCheckCommandLine(option)

  realization%waypoint_list => WaypointListCreate()
  
  ! initialize flow mode
  if (len_trim(option%flowmode) > 0) then
    ! set the operational mode (e.g.  MPH_MODE, etc)
    call InitSetFlowMode(option)
  else
    option%nphase = 1
    option%liquid_phase = 1
    option%use_isothermal = PETSC_TRUE  ! assume default isothermal when only transport
    call TimestepperDestroy(simulation%flow_timestepper)
  endif
    
  ! initialize transport mode
  if (option%ntrandof > 0) then
  else
    call TimestepperDestroy(simulation%tran_timestepper)
  endif

  ! initialize surface-flow mode
  if (option%surf_flow_on) then
    call setSurfaceFlowMode(option)
    surf_realization%waypoint_list => WaypointListCreate()
  else
    call TimestepperDestroy(simulation%surf_flow_timestepper)
  endif

  ! initialize surface-flow mode
  if (option%ngeomechdof > 0) then
    geomech_realization%waypoint_list => WaypointListCreate()
  else
    call TimestepperDestroy(simulation%geomech_timestepper)
  endif

  ! read in the remainder of the input file
  call InitReadInput(simulation)
  call InputDestroy(realization%input)

  ! destroy other 'input' used above
  ! (TODO: fmy- need to check if this is the right place to do so)
  call InputDestroy(surf_realization%input)    !allocated in Line 182)
  call InputDestroy(geomech_realization%input) ! allocated in Line 186)

end subroutine Init

! ************************************************************************** !

subroutine InitReadInputFilenames(option,filenames)
  ! 
  ! Reads filenames for multi-simulation runs
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 

  use Option_module
  use Input_Aux_module

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: filename_count
  type(input_type), pointer :: input
  PetscBool :: card_found

  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = "FILENAMES"
  call InputFindStringInFile(input,option,string) 

  card_found = PETSC_FALSE
  if (InputError(input)) then
    ! if the FILENAMES card is not included, we will assume that only
    ! filenames exist in the file.
    rewind(input%fid)
  else
    card_found = PETSC_TRUE
  endif
    
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
  enddo
  
  allocate(filenames(filename_count))
  filenames = ''
  rewind(input%fid) 

  if (card_found) then
    string = "FILENAMES"
    call InputFindStringInFile(input,option,string) 
  endif
  
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
    filenames(filename_count) = filename
  enddo

  call InputDestroy(input)

end subroutine InitReadInputFilenames

! ************************************************************************** !

subroutine InitSubsurfaceReadRequiredCards(realization)
  ! 
  ! Reads pflow input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07, refactored 08/20/14
  ! 

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_class
  use HDF5_Aux_module

  use General_module
  use Reaction_module  
  use Reaction_Aux_module  

  implicit none

  type(realization_type) :: realization

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch, patch2 
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization
  
  input => realization%input
  
! Read in select required cards
!.........................................................................
  
  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call DiscretizationReadRequiredCards(discretization,input,option)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(realization%patch_list)) then
        realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,realization%patch_list)
      realization%patch => patch
  end select

  ! optional required cards - yes, an oxymoron, but we need to know if
  ! these exist before we can go any further.
  rewind(input%fid)  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    select case(trim(card))

!....................
#ifndef INIT_REFACTOR    
      case ('MODE')
        call InputReadWord(input,option,option%flowmode,PETSC_TRUE)
        call InputErrorMsg(input,option,'flowmode','mode')
        select case(trim(option%flowmode))
          case('GENERAL')
            call GeneralRead(input,option)
        end select
#endif  
!....................
      case('DBASE_FILENAME')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'filename','DBASE_FILENAME')
        if (index(word,'.h5') > 0) then
#if defined(PETSC_HAVE_HDF5)
          call HDF5ReadDbase(word,option)
#endif
        else
          call InputReadASCIIDbase(word,option)
        endif

!....................
#if defined(SCORPIO)
      case('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
        call InputSkipToEnd(input,option,'HDF5_WRITE_GROUP_SIZE')

      case('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')
#endif

!....................
      case('PROC')
        ! processor decomposition
        if (realization%discretization%itype == STRUCTURED_GRID) then
          grid => realization%patch%grid
          ! strip card from front of string
          call InputReadInt(input,option,grid%structured_grid%npx)
          call InputDefaultMsg(input,option,'npx')
          call InputReadInt(input,option,grid%structured_grid%npy)
          call InputDefaultMsg(input,option,'npy')
          call InputReadInt(input,option,grid%structured_grid%npz)
          call InputDefaultMsg(input,option,'npz')
 
          if (option%myrank == option%io_rank .and. &
              option%print_to_screen) then
            option%io_buffer = ' Processor Decomposition:'
            call printMsg(option)
            write(option%io_buffer,'("  npx   = ",3x,i4)') &
              grid%structured_grid%npx
            call printMsg(option)
            write(option%io_buffer,'("  npy   = ",3x,i4)') &
              grid%structured_grid%npy
            call printMsg(option)
            write(option%io_buffer,'("  npz   = ",3x,i4)') &
              grid%structured_grid%npz
            call printMsg(option)
          endif
  
          if (option%mycommsize /= grid%structured_grid%npx * &
                                 grid%structured_grid%npy * &
                                 grid%structured_grid%npz) then
            write(option%io_buffer,*) 'Incorrect number of processors specified: ', &
                           grid%structured_grid%npx*grid%structured_grid%npy* &
                           grid%structured_grid%npz,' commsize = ',option%mycommsize
            call printErrMsg(option)
          endif
        endif
  
!....................
      case('CHEMISTRY')
        !geh: for some reason, we need this with CHEMISTRY read for 
        !     multicontinuum
        option%use_mc = PETSC_TRUE
        call ReactionInit(realization%reaction,input,option)
    end select
  enddo
  
#if defined(SCORPIO)
  call InitCommonCreateIOGroups(option)
#endif  

end subroutine InitSubsurfaceReadRequiredCards

! ************************************************************************** !

subroutine InitReadInput(simulation)
  ! 
  ! Reads pflow input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Simulation_module
  use Option_module
  use Field_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Structured_module
  use Solver_module
  use Material_module
  use Saturation_Function_module  
  use Characteristic_Curves_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Fluid_module
  use Realization_class
  use Timestepper_module
  use Region_module
  use Condition_module
  use Transport_Constraint_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Integral_Flux_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Reaction_Aux_module
  use Discretization_module
  use Input_Aux_module
  use String_module
  use Units_module
  use Uniform_Velocity_module
  use Reaction_Mineral_module
  use Regression_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Mass_Transfer_module
  use EOS_module
  use EOS_Water_module
  use SrcSink_Sandbox_module
  use Creep_Closure_module
  
  use Surface_Flow_module
  use Surface_Init_module, only : SurfaceInitReadInput
  use Geomechanics_Init_module, only : GeomechanicsInitReadInput
  use Geomechanics_Realization_class
#ifdef SOLID_SOLUTION
  use Reaction_Solid_Solution_module, only : SolidSolutionReadFromInputFile
#endif
 
  implicit none
  
  type(simulation_type) :: simulation

  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
    
  PetscBool :: continuation_flag
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal :: units_conversion
  PetscInt :: temp_int
  PetscInt :: count, id
  
  PetscBool :: vel_cent
  PetscBool :: vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate
  
  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_type), pointer :: sec_tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  type(integral_flux_type), pointer :: integral_flux
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_property_type), pointer :: material_property
  type(fluid_property_type), pointer :: fluid_property
  type(saturation_function_type), pointer :: saturation_function
  class(characteristic_curves_type), pointer :: characteristic_curves

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch   
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  type(solver_type), pointer :: default_solver
  type(solver_type), pointer :: solver_pointer
  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
  type(timestepper_type), pointer :: default_timestepper
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(uniform_velocity_dataset_type), pointer :: uniform_velocity_dataset
  class(dataset_base_type), pointer :: dataset
  type(mass_transfer_type), pointer :: flow_mass_transfer
  type(mass_transfer_type), pointer :: rt_mass_transfer
  type(input_type), pointer :: input
  type(geomech_realization_type), pointer :: geomech_realization

  nullify(flow_timestepper)
  nullify(tran_timestepper)
  nullify(flow_solver)
  nullify(tran_solver)
  
  realization => simulation%realization
  patch => realization%patch
  
  geomech_realization => simulation%geomech_realization

  if (associated(patch)) grid => patch%grid

  option => realization%option
  output_option => realization%output_option
  field => realization%field
  reaction => realization%reaction
  input => realization%input

  tran_timestepper => simulation%tran_timestepper
  if (associated(tran_timestepper)) then
    tran_solver => tran_timestepper%solver
    tran_solver%itype = TRANSPORT_CLASS
  endif
  flow_timestepper => simulation%flow_timestepper
  if (associated(flow_timestepper)) then
    flow_solver => flow_timestepper%solver
    flow_solver%itype = FLOW_CLASS
  endif

  if (associated(flow_timestepper)) then
    default_timestepper => flow_timestepper
    default_solver => flow_solver
  else
    default_timestepper => tran_timestepper
    default_solver => tran_solver
  endif

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                              
  rewind(input%fid)  
      
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))

!....................
      case ('MODE')
         call InputReadWord(input, option, word, PETSC_FALSE)
         call StringToUpper(word)
         if ('TH' == trim(word)) then
            call InputReadWord(input, option, word, PETSC_TRUE)
            call InputErrorMsg(input, option, 'th freezing mode', 'mode th')
            call StringToUpper(word)
            if ('FREEZING' == trim(word)) then
               option%use_th_freezing = PETSC_TRUE
               option%io_buffer = ' TH: using FREEZING submode!'
               call printMsg(option)
               ! Override the default setting for TH-mode with freezing
               call EOSWaterSetDensityPainter()
               call EOSWaterSetEnthalpyPainter()
            else if ('NO_FREEZING' == trim(word)) then
               option%use_th_freezing = PETSC_FALSE
               option%io_buffer = ' TH: using NO_FREEZING submode!'
               call printMsg(option)
            else
               ! NOTE(bja, 2013-12) use_th_freezing defaults to false, can skip this....
               option%io_buffer = ' TH: must specify FREEZING or NO_FREEZING submode!'
               call printErrMsg(option)
            endif
         else if (trim(word) == 'GENERAL') then
           call InputReadWord(input, option, word, PETSC_TRUE)
           if (input%ierr == 0) then
             call InputSkipToEnd(input,option,card)
           endif
         endif  
         
!....................
      case('CREEP_CLOSURE')
        call CreepClosureInit()
        creep_closure => CreepClosureCreate()
        call creep_closure%Read(input,option)
        option%flow%transient_porosity = PETSC_TRUE
        
!....................
      case ('ICE_MODEL')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('PAINTER_EXPLICIT')
            option%ice_model = PAINTER_EXPLICIT
          case ('PAINTER_KARRA_IMPLICIT')
            option%ice_model = PAINTER_KARRA_IMPLICIT
          case ('PAINTER_KARRA_EXPLICIT')
            option%ice_model = PAINTER_KARRA_EXPLICIT
          case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
            option%ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
          case ('DALL_AMICO')
            option%ice_model = DALL_AMICO
          case default
            option%io_buffer = 'Cannot identify the specificed ice model.' // &
             'Specify PAINTER_EXPLICIT or PAINTER_KARRA_IMPLICIT' // &
             ' or PAINTER_KARRA_EXPLICIT or PAINTER_KARRA_EXPLICIT_NOCRYO ' // &
             ' or DALL_AMICO.'
            call printErrMsg(option)
          end select

!....................
      case ('ONLY_VERTICAL_FLOW')
        option%flow%only_vertical_flow = PETSC_TRUE
        if (option%iflowmode /= TH_MODE) then
          option%io_buffer = 'ONLY_VERTICAL_FLOW implemented in TH_MODE'
          call printErrMsg(option)
        endif

!....................
      case ('RELATIVE_PERMEABILITY_AVERAGE')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('UPWIND')
            option%rel_perm_aveg = UPWIND
          case ('HARMONIC')
            option%rel_perm_aveg = HARMONIC
          case ('DYNAMIC_HARMONIC')
            option%rel_perm_aveg = DYNAMIC_HARMONIC
          case default
            option%io_buffer = 'Cannot identify the specificed ' // &
              'RELATIVE_PERMEABILITY_AVERAGE.'
            call printErrMsg(option)
          end select

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('CHEMISTRY')
        call ReactionReadPass2(reaction,input,option)

!....................
      case ('NONUNIFORM_VELOCITY')
        call InputReadNChars(input,option, &
                             realization%nonuniform_velocity_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','NONUNIFORM_VELOCITY')

      case ('UNIFORM_VELOCITY')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        uniform_velocity_dataset%rank = 3
        uniform_velocity_dataset%interpolation_method = 1 ! 1 = STEP
        uniform_velocity_dataset%is_cyclic = PETSC_FALSE
        allocate(uniform_velocity_dataset%times(1))
        uniform_velocity_dataset%times = 0.d0
        allocate(uniform_velocity_dataset%values(3,1))
        uniform_velocity_dataset%values = 0.d0
        call InputReadDouble(input,option,uniform_velocity_dataset%values(1,1))
        call InputErrorMsg(input,option,'velx','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(2,1))
        call InputErrorMsg(input,option,'vely','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(3,1))
        call InputErrorMsg(input,option,'velz','UNIFORM_VELOCITY')
        ! read units, if present
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          units_conversion = UnitsConvertToInternal(word,option) 
          uniform_velocity_dataset%values(:,1) = &
            uniform_velocity_dataset%values(:,1) * units_conversion
        endif
        call UniformVelocityDatasetVerify(option,uniform_velocity_dataset)
        realization%uniform_velocity_dataset => uniform_velocity_dataset
      
      case ('VELOCITY_DATASET')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        call UniformVelocityDatasetRead(uniform_velocity_dataset,input,option)
        realization%uniform_velocity_dataset => uniform_velocity_dataset

!....................
      case ('DEBUG')
        call DebugRead(realization%debug,input,option)
               
!....................
      case ('PRINT_PRIMAL_GRID')
        option%print_explicit_primal_grid = PETSC_TRUE
        
!....................
      case ('PRINT_DUAL_GRID')
        option%print_explicit_dual_grid = PETSC_TRUE

!....................
      case ('MAX_CHANGE')
        call InputReadDouble(input,option,option%dpmxe)
        call InputErrorMsg(input,option,'dpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dtmpmxe)
        call InputErrorMsg(input,option,'dtmpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dsmxe)
        call InputErrorMsg(input,option,'dsmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dcmxe)
        call InputErrorMsg(input,option,'dcmxe','MAX_CHANGE')
        
!....................
      case ('PROC')
      
!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION') 
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%region_list)   
        nullify(region)   

!....................
      case ('FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'FLOW_CONDITION','name') 
        call printMsg(option,flow_condition%name)
        if (option%iflowmode == G_MODE) then
          call FlowConditionGeneralRead(flow_condition,input,option)
        else
          call FlowConditionRead(flow_condition,input,option)
        endif
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)
        
!....................
      case ('TRANSPORT_CONDITION')
        if (.not.associated(reaction)) then
          option%io_buffer = 'TRANSPORT_CONDITIONs not supported without ' // &
            'CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_condition => TranConditionCreate(option)
        call InputReadWord(input,option,tran_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'TRANSPORT_CONDITION','name') 
        call printMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition,realization%transport_constraints, &
                               reaction,input,option)
        call TranConditionAddToList(tran_condition,realization%transport_conditions)
        nullify(tran_condition)

!....................
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name') 
        call printMsg(option,tran_constraint%name)
        call TranConstraintRead(tran_constraint,reaction,input,option)
        call TranConstraintAddToList(tran_constraint,realization%transport_constraints)
        nullify(tran_constraint)


!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('SOURCE_SINK_SANDBOX')
        call SSSandboxInit(option)
        call SSSandboxRead(input,option)
      
!....................
      case ('FLOW_MASS_TRANSFER')
        flow_mass_transfer => MassTransferCreate()
        call InputReadWord(input,option,flow_mass_transfer%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Flow Mass Transfer name') 
        call MassTransferRead(flow_mass_transfer,input,option)
        call MassTransferAddToList(flow_mass_transfer, &
                                   realization%flow_mass_transfer_list)
        nullify(flow_mass_transfer)
      
!....................
      case ('RT_MASS_TRANSFER')
        rt_mass_transfer => MassTransferCreate()
        call InputReadWord(input,option,rt_mass_transfer%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'RT Mass Transfer name')
        call MassTransferRead(rt_mass_transfer,input,option)
        call MassTransferAddToList(rt_mass_transfer, &
                                   realization%rt_mass_transfer_list)
        nullify(rt_mass_transfer)
      
!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization,strata)
        nullify(strata)

!.....................
      case ('DATASET')
        nullify(dataset)
        call DatasetRead(input,dataset,option)
        call DatasetBaseAddToList(dataset,realization%datasets)
        nullify(dataset)
        
!....................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_pressure)
        call InputDefaultMsg(input,option,'Reference Pressure') 

!....................

      case('REFERENCE_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_water_density)
        call InputDefaultMsg(input,option,'Reference Density') 

!....................

      case('MINIMUM_HYDROSTATIC_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%minimum_hydrostatic_pressure)
        call InputDefaultMsg(input,option,'Minimum Hydrostatic Pressure') 

!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_temperature)
        call InputDefaultMsg(input,option,'Reference Temperature') 

!......................

      case('REFERENCE_POROSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_porosity)
        call InputDefaultMsg(input,option,'Reference Porosity') 

!......................

      case('REFERENCE_SATURATION')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_saturation)
        call InputDefaultMsg(input,option,'Reference Saturation') 

!......................

      case('NONISOTHERMAL')
        option%use_isothermal = PETSC_FALSE

!......................

      case('ISOTHERMAL')
        option%use_isothermal = PETSC_TRUE
        
!......................

      case('UPDATE_FLOW_PERMEABILITY')
        option%update_flow_perm = PETSC_TRUE
        
!......................

      case('DFN')
        grid%unstructured_grid%grid_type = TWO_DIM_GRID    
            
!......................

      case("MULTIPLE_CONTINUUM")
        option%use_mc = PETSC_TRUE
              
!......................

      case('SECONDARY_CONTINUUM_SOLVER')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can only be used ' // &
                             'with MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif      
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('KEARST')
            option%secondary_continuum_solver = 1
          case('HINDMARSH')
            option%secondary_continuum_solver = 2
          case('THOMAS')
            option%secondary_continuum_solver = 3
          case default
            option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                               'HINDMARSH or KEARST. For single component'// &
                               'chemistry THOMAS can be used.'
          call printErrMsg(option)    
        end select        
!....................

      case('SECONDARY_CONSTRAINT')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONSTRAINT can only be used with ' // &
                             'MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif
        if (.not.associated(reaction)) then
          option%io_buffer = 'SECONDARY_CONSTRAINT not supported without' // &
                             'CHEMISTRY.'
          call printErrMsg(option)
        endif
        sec_tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,sec_tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'secondary constraint','name') 
        call printMsg(option,sec_tran_constraint%name)
        call TranConstraintRead(sec_tran_constraint,reaction,input,option)
        realization%sec_transport_constraint => sec_tran_constraint
        nullify(sec_tran_constraint)        

!......................

      case('BRIN','BRINE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%m_nacl)
        call InputDefaultMsg(input,option,'NaCl Concentration') 

        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            !strictly correct if no CO2 present initially-pcl
            option%m_nacl = option%m_nacl /FMWNACL/(1.D0-option%m_nacl)
          case('MOLE')    
            !strictly correct if no CO2 present initially-pcl
            option%m_nacl = option%m_nacl /FMWH2O/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (OptionPrintToScreen(option)) print *, option%m_nacl
!......................

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call InputReadNChars(input,option,option%restart_filename,MAXSTRINGLENGTH, &
                             PETSC_TRUE)
        call InputErrorMsg(input,option,'RESTART','Restart file name') 
        call InputReadDouble(input,option,option%restart_time)
        if (input%ierr == 0) then
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) then
            option%restart_time = option%restart_time* &
                                  UnitsConvertToInternal(word,option)
          else
            call InputDefaultMsg(input,option,'RESTART, time units')
          endif
        endif

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call InputReadInt(input,option,option%checkpoint_frequency)

        if (input%ierr == 1) then
          option%checkpoint_frequency = 0
          do
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option,card)
            if (InputCheckExit(input,option)) exit

            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','CHECKPOINT')
            call StringToUpper(word)

            select case(trim(word))
              case ('PERIODIC')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'time increment', &
                                   'OUTPUT,PERIODIC')
                call StringToUpper(word)

                select case(trim(word))
                  case('TIME')
                    call InputReadDouble(input,option,temp_real)
                    call InputErrorMsg(input,option,'time increment', &
                                       'CHECKPOINT,PERIODIC,TIME')
                    call InputReadWord(input,option,word,PETSC_TRUE)
                    call InputErrorMsg(input,option,'time increment units', &
                                       'CHECKPOINT,PERIODIC,TIME')
                    units_conversion = UnitsConvertToInternal(word,option)
                    output_option%periodic_checkpoint_time_incr = temp_real* &
                                                              units_conversion
                  case('TIMESTEP')
                    call InputReadInt(input,option,option%checkpoint_frequency)
                    call InputErrorMsg(input,option,'timestep increment', &
                                       'CHECKPOINT,PERIODIC,TIMESTEP')
                  case default
                    option%io_buffer = 'Keyword: ' // trim(word) // &
                                       ' not recognized in CHECKPOINT,PERIODIC.'
                    call printErrMsg(option)
                end select
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                                   ' not recognized in CHECKPOINT.'
                call printErrMsg(option)
            end select
          enddo
          if (output_option%periodic_checkpoint_time_incr /= 0.d0 .and. &
              option%checkpoint_frequency /= 0) then
            option%io_buffer = 'Both TIME and TIMESTEP cannot be specified ' // &
              'for CHECKPOINT,PERIODIC.'
            call printErrMsg(option)
          endif
          if (output_option%periodic_checkpoint_time_incr == 0.d0 .and. &
              option%checkpoint_frequency == 0) then
            option%io_buffer = 'Either, TIME and TIMESTEP need to be specified ' // &
              'for CHECKPOINT,PERIODIC.'
            call printErrMsg(option)
          endif
        endif

!......................

      case ('NUMERICAL_JACOBIAN_FLOW')
        option%numerical_derivatives_flow = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_RXN')
        option%numerical_derivatives_rxn = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_MULTI_COUPLE')
        option%numerical_derivatives_multi_coupling = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS')
        option%compute_statistics = PETSC_TRUE

!....................

      case ('CO2_DATABASE')
        call InputReadNChars(input,option,option%co2_database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'CO2_DATABASE','filename')
        
!....................

      case ('TIMESTEPPER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            if (associated(flow_solver)) then
              call TimestepperRead(flow_timestepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call TimestepperRead(tran_timestepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case default
            option%io_buffer = 'TIMESTEPPER must specify FLOW or TRANSPORT.'
            call printErrMsg(option)
            if (associated(default_timestepper)) then
              call TimestepperRead(default_timestepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
        end select

!....................

      case ('LINEAR_SOLVER')
        nullify(solver_pointer)
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            solver_pointer => flow_solver
          case('TRAN','TRANSPORT')
            solver_pointer => tran_solver
          case default
            solver_pointer => default_solver
        end select
        if (associated(solver_pointer)) then
          call SolverReadLinear(solver_pointer,input,option)
        else
          call InputSkipToEnd(input,option,card)
        endif

!....................

      case ('NEWTON_SOLVER')
        nullify(solver_pointer)
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            solver_pointer => flow_solver
          case('TRAN','TRANSPORT')
            solver_pointer => tran_solver
          case default
            solver_pointer => default_solver
        end select
        if (associated(solver_pointer)) then
          call SolverReadNewton(solver_pointer,input,option)
        else
          call InputSkipToEnd(input,option,card)
        endif

        if (associated(solver_pointer,flow_solver) .and. &
            solver_pointer%check_post_convergence) then
          option%flow%check_post_convergence = PETSC_TRUE
          option%flow%inf_scaled_res_tol = &
            solver_pointer%newton_inf_scaled_res_tol
          option%flow%inf_rel_update_tol = &
            solver_pointer%newton_inf_rel_update_tol
        endif
        if (associated(solver_pointer,tran_solver) .and. &
            solver_pointer%check_post_convergence) then
          option%transport%check_post_convergence = PETSC_TRUE
          option%transport%inf_scaled_res_tol = &
            solver_pointer%newton_inf_scaled_res_tol
          option%transport%inf_rel_update_tol = &
            solver_pointer%newton_inf_rel_update_tol
        endif
        
!....................

      case ('FLUID_PROPERTY')

        fluid_property => FluidPropertyCreate()
        call FluidPropertyRead(fluid_property,input,option)
        call FluidPropertyAddToList(fluid_property,realization%fluid_properties)
        nullify(fluid_property)
        
!....................

      case ('EOS')
        call EOSRead(input,option)

!....................

      case ('SATURATION_FUNCTION')
        saturation_function => SaturationFunctionCreate(option)
        call InputReadWord(input,option,saturation_function%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SATURATION_FUNCTION')
        call SaturationFunctionRead(saturation_function,input,option)
        call SatFunctionComputePolynomial(option,saturation_function)
        call PermFunctionComputePolynomial(option,saturation_function)
        call SaturationFunctionAddToList(saturation_function, &
                                         realization%saturation_functions)
        nullify(saturation_function)   

!....................

      case ('CHARACTERISTIC_CURVES')
      
        characteristic_curves => CharacteristicCurvesCreate()
        call InputReadWord(input,option,characteristic_curves%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','CHARACTERISTIC_CURVES')
        call CharacteristicCurvesRead(characteristic_curves,input,option)
!        call SatFunctionComputePolynomial(option,saturation_function)
!        call PermFunctionComputePolynomial(option,saturation_function)
        call CharacteristicCurvesAddToList(characteristic_curves, &
                                          realization%characteristic_curves)
        nullify(characteristic_curves)   

!....................
      
      case ('MATERIAL_PROPERTY')

        material_property => MaterialPropertyCreate()
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')        
        call MaterialPropertyRead(material_property,input,option)
        call MaterialPropertyAddToList(material_property,realization%material_properties)
        nullify(material_property)

!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
!        call PetscOptionsInsertString('-viewer_binary_mpiio')

      case ('HANDSHAKE_IO')
        call InputReadInt(input,option,option%io_handshake_buffer_size)
        call InputErrorMsg(input,option,'io_handshake_buffer_size','HANDSHAKE_IO')

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%overwrite_restart_transport = PETSC_TRUE

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%overwrite_restart_flow = PETSC_TRUE

      case ('INITIALIZE_FLOW_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_flow_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_FLOW_FROM_FILE') 

      case ('INITIALIZE_TRANSPORT_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_transport_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_TRANSPORT_FROM_FILE') 

      case ('CENTRAL_DIFFERENCE')
        option%use_upwinding = PETSC_FALSE

!....................
      case ('OBSERVATION')
        observation => ObservationCreate()
        call ObservationRead(observation,input,option)
        call ObservationAddToList(observation, &
                                  realization%patch%observation_list)
        nullify(observation)
      
!....................
      case ('INTEGRAL_FLUX')
        integral_flux => IntegralFluxCreate()
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Integral Flux name') 
        call IntegralFluxRead(integral_flux,input,option)
        call IntegralFluxAddToList(integral_flux, &
                                   realization%patch%integral_flux_list)
        nullify(integral_flux)
      
!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time','WALLCLOCK_STOP') 

        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) word = 'h'
        call InputDefaultMsg(input,option,'WALLCLOCK_STOP time units')
        units_conversion = UnitsConvertToInternal(word,option) 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*units_conversion
      
!....................
      case ('OUTPUT')
        vel_cent = PETSC_FALSE
        vel_face = PETSC_FALSE
        fluxes = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','OUTPUT') 
          call StringToUpper(word)
          select case(trim(word))
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial = PETSC_FALSE
            case('PROCESSOR_ID')
              option%io_buffer = 'PROCESSOR_ID output must now be entered under OUTPUT/VARIABLES card as PROCESS_ID.'
              call printErrMsg(option)
!              output_option%print_iproc = PETSC_TRUE
            case('PERMEABILITY')
              option%io_buffer = 'PERMEABILITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              option%io_buffer = 'POROSITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_porosity = PETSC_TRUE
            case('TORTUOSITY')
              option%io_buffer = 'TORTUOSITY output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_tortuosity = PETSC_TRUE
            case('VOLUME')
              option%io_buffer = 'VOLUME output must now be entered under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_volume = PETSC_TRUE
            case('MASS_BALANCE')
              option%compute_mass_balance_new = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputDefaultMsg(input,option, &
                                 'MASS_BALANCE,DETAILED,OUTPUT')
              if (len_trim(word) > 0) then
                call StringToUpper(word)
                select case(trim(word))
                  case('DETAILED')
                    option%mass_bal_detailed = PETSC_TRUE
                  case('DEFAULT')
                    option%io_buffer = 'Keyword: ' // trim(word) // &
                      ' not recognized in OUTPUT,'// &
                      'MASS_BALANCE,DETAILED.'
                    call printErrMsg(option)
                end select
              endif
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE
            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','OUTPUT')
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
                    call WaypointInsertInList(waypoint,realization%waypoint_list)
                  endif
                enddo
                if (.not.continuation_flag) exit
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
              enddo
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,OUTPUT_FILE')
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
                                     'OUTPUT,PERIODIC,SCREEN')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,SCREEN.'
                  call printErrMsg(option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word,option) 
                  output_option%periodic_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then

                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word,option) 
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_output = PETSC_TRUE    
                        call WaypointInsertInList(waypoint,realization%waypoint_list)
                        temp_real = temp_real + output_option%periodic_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'OUTPUT,PERIODIC,TIME')
                    endif
                  endif                  
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,TIMESTEP')
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
                  if (len_trim(word) > 0) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadWord(input,option,word,PETSC_TRUE)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option,'timestep increment', &
                                        'OUTPUT,FORMAT,MULTIPLE_FILES,TIMES_PER_FILE')
                            case default
                              option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,'// &
                                     'FORMAT,MULTIPLE_FILES,TIMES_PER_FILE.'
                              call printErrMsg(option)
                          end select
                        endif
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recognized.  Use "SINGLE_FILE" or ' // &
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
                    case('FEBRICK')
                      output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
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
                    output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,FORMAT.'
                  call printErrMsg(option)
              end select
            case('VELOCITY_AT_CENTER')
              vel_cent = PETSC_TRUE
            case('VELOCITY_AT_FACE')
              vel_face = PETSC_TRUE
            case('FLUXES')
              fluxes = PETSC_TRUE
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
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
            case('VARIABLES')
              call OutputVariableRead(input,option,output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputVariableRead(input,option,output_option%aveg_output_variable_list)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in OUTPUT.'
              call printErrMsg(option)              
          end select
        enddo
        if (vel_cent) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_cent = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_vel_cent = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_vel_cent = PETSC_TRUE
        endif
        if (vel_face) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_face = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_vel_face = PETSC_TRUE
        endif
        if (fluxes) then
          output_option%print_fluxes = PETSC_TRUE
        endif
        if(output_option%aveg_output_variable_list%nvars>0) then
          if(output_option%periodic_output_time_incr==0.d0) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without' // &
                               ' PERIODIC TIME being set.'
            call printErrMsg(option)
          endif
          if(.not.output_option%print_hdf5) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for FORMAT HDF5'
            call printErrMsg(option)
          endif
        endif
        if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate.or.aveg_energy_flowrate) then
          if (output_option%print_hdf5) then
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
           option%flow%store_fluxes = PETSC_TRUE
          endif
          if (associated(grid%unstructured_grid%explicit_grid)) then
           option%flow%store_fluxes = PETSC_TRUE
            output_option%print_explicit_flowrate = mass_flowrate
          endif
        
        endif

!.....................
      case ('REGRESSION')
        call RegressionRead(simulation%regression,input,option)

!.....................
      case ('TIME')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','TIME') 
          select case(trim(word))
            case('STEADY_STATE')
              option%steady_state = PETSC_TRUE
            case('FINAL_TIME')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units','TIME')
              realization%output_option%tunit = trim(word)
              realization%output_option%tconv = UnitsConvertToInternal(word,option)
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*realization%output_option%tconv
              waypoint%print_output = PETSC_TRUE              
              call WaypointInsertInList(waypoint,realization%waypoint_list)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Initial Timestep Size Time Units','TIME')
              default_timestepper%dt_init = temp_real*UnitsConvertToInternal(word,option)
            case('MINIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Minimum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Minimum Timestep Size Time Units','TIME')
              default_timestepper%dt_min = temp_real*UnitsConvertToInternal(word,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Maximum Timestep Size Time Units','TIME')
              waypoint => WaypointCreate()
              waypoint%dt_max = temp_real*UnitsConvertToInternal(word,option)
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time','TIME') 
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time Units','TIME')
                  waypoint%time = temp_real*UnitsConvertToInternal(word,option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" after ' // &
                                     'maximum timestep size should be "at".'
                  call printErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif     
              call WaypointInsertInList(waypoint,realization%waypoint_list)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in TIME.'
              call printErrMsg(option)              
          end select
        enddo

        if (associated(flow_timestepper)) then
          flow_timestepper%dt_init = default_timestepper%dt_init
        endif
        if (associated(tran_timestepper)) then
          tran_timestepper%dt_init = default_timestepper%dt_init
        endif
        option%flow_dt = default_timestepper%dt_init
        option%tran_dt = default_timestepper%dt_init
      
!.....................
      case ('SURFACE_FLOW')
        call SurfaceInitReadInput(simulation%surf_realization, &
                              simulation%surf_flow_timestepper%solver,input,option)
        simulation%surf_flow_timestepper%dt_init = simulation%surf_realization%dt_init
        simulation%surf_flow_timestepper%dt_max = simulation%surf_realization%dt_max
        option%surf_subsurf_coupling_flow_dt = simulation%surf_realization%dt_coupling
        option%surf_flow_dt=simulation%surf_flow_timestepper%dt_init

        ! Add first waypoint
        waypoint => WaypointCreate()
        waypoint%time = 0.d0
        call WaypointInsertInList(waypoint,simulation%surf_realization%waypoint_list)

        ! Add final_time waypoint to surface_realization
        waypoint => WaypointCreate()
        waypoint%final = PETSC_TRUE
        waypoint%time = realization%waypoint_list%last%time
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,simulation%surf_realization%waypoint_list)

!......................
      case ('GEOMECHANICS')
        call GeomechanicsInitReadInput(geomech_realization, &
                         simulation%geomech_timestepper%solver,input,option)
        ! Add first waypoint
        waypoint => WaypointCreate()
        waypoint%time = 0.d0
        call WaypointInsertInList(waypoint,simulation%geomech_realization%waypoint_list)

        ! Add final_time waypoint to geomech_realization
        waypoint => WaypointCreate()
        waypoint%final = PETSC_TRUE
        waypoint%time = realization%waypoint_list%last%time
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,simulation%geomech_realization%waypoint_list)

!......................
      case ('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!......................
      case ('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')

!....................
      case ('DBASE_FILENAME')

!....................
      case default
    
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select

  enddo
                                      
end subroutine InitReadInput

! ************************************************************************** !

subroutine InitSetFlowMode(option)
  ! 
  ! Sets the flow mode (richards, vadose, mph, etc.)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  use Option_module
  use String_module
  use General_Aux_module

  implicit none 

  type(option_type) :: option
  
  call StringToUpper(option%flowmode)
  select case(option%flowmode)
    case('TH')
      option%iflowmode = TH_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
      option%flow%store_fluxes = PETSC_TRUE
    case('MIS','MISCIBLE')
      option%iflowmode = MIS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 2
      option%nflowspec = 2
      option%io_buffer = 'Material Auxvars must be refactored for MISCIBLE.'
      call printErrMsg(option)
    case('RICHARDS')
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    case('MPH','MPHASE')
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2 ! read CO2DATA0.dat
!     option%itable = 1 ! create CO2 database: co2data.dat
      option%use_isothermal = PETSC_FALSE
    case('FLASH2')
      option%iflowmode = FLASH2_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%use_isothermal = PETSC_FALSE
   case('IMS','IMMIS','THS')
      option%iflowmode = IMS_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%io_buffer = 'Material Auxvars must be refactored for IMMIS.'
      call printErrMsg(option)
    case('GENERAL')
      option%iflowmode = G_MODE
      option%nphase = 2
      option%liquid_phase = 1  ! liquid_pressure
      option%gas_phase = 2     ! gas_pressure

      option%air_pressure_id = 3
      option%capillary_pressure_id = 4
      option%vapor_pressure_id = 5
      option%saturation_pressure_id = 6

      option%water_id = 1
      option%air_id = 2
      option%energy_id = 3

      option%nflowdof = 3
      option%nflowspec = 2
      option%use_isothermal = PETSC_FALSE
    case default
      option%io_buffer = 'Mode: '//trim(option%flowmode)//' not recognized.'
      call printErrMsg(option)
  end select
  
end subroutine InitSetFlowMode

! ************************************************************************** !

subroutine setSurfaceFlowMode(option)
  ! 
  ! Sets the flow mode for surface (richards, th, etc.)
  ! 
  ! Author: Gautam Bisht
  ! Date: 07/30/14
  ! 

  use Option_module
  use String_module

  implicit none 

  type(option_type) :: option
  
  call StringToUpper(option%flowmode)
  select case(option%flowmode)
    case('RICHARDS')
      option%nsurfflowdof = ONE_INTEGER
    case('TH')
      option%nsurfflowdof = TWO_INTEGER
    case default
      option%io_buffer = 'Mode: '//trim(option%flowmode)//' not recognized.'
      call printErrMsg(option)
  end select
  
end subroutine setSurfaceFlowMode

! ************************************************************************** !

subroutine InitCommonVerifyAllCouplers(realization)
  ! 
  ! Verifies the connectivity of a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/8/08
  ! 

  use Realization_class
  use Patch_module
  use Coupler_module

  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%initial_condition_list)
      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%boundary_condition_list)
      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%source_sink_list)

    cur_patch => cur_patch%next
  enddo
  
end subroutine InitCommonVerifyAllCouplers

! ************************************************************************** !

subroutine InitCommonVerifyCoupler(realization,patch,coupler_list)
  ! 
  ! Verifies the connectivity of a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/8/08
  ! 

  use Realization_class
  use Discretization_module
  use Option_module 
  use Coupler_module
  use Condition_module
  use Grid_module
  use Output_module
  use Output_Tecplot_module, only : OutputVectorTecplot  
  use Patch_module

  implicit none

  type(realization_type) :: realization
  type(coupler_list_type), pointer :: coupler_list

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(coupler_type), pointer :: coupler
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: iconn, icell, local_id
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (.not.associated(coupler_list)) return

  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
    if (associated(coupler%connection_set)) then
      do iconn = 1, coupler%connection_set%num_connections
        local_id = coupler%connection_set%id_dn(iconn)
!        vec_ptr(local_id) = coupler%id
!geh: let's sum the # of connections
         vec_ptr(local_id) = vec_ptr(local_id) + 1
      enddo
    else
      if (associated(coupler%region)) then
        do icell = 1, coupler%region%num_cells
          local_id = coupler%region%cell_ids(icell)
!          vec_ptr(local_id) = coupler%id
         vec_ptr(local_id) = vec_ptr(local_id) + 1
        enddo
      endif
    endif
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
    if (len_trim(coupler%flow_condition_name) > 0) then
      dataset_name = coupler%flow_condition_name
    else if (len_trim(coupler%tran_condition_name) > 0) then
      dataset_name = coupler%tran_condition_name
    endif
    write(word,*) patch%id
    dataset_name = trim(dataset_name) // '_' // &
                   trim(coupler%region%name) // '_' // &
                   trim(adjustl(word))
    dataset_name = dataset_name(1:28)
    filename = trim(dataset_name) // '.tec'
    call OutputVectorTecplot(filename,dataset_name,realization,global_vec)

    coupler => coupler%next
  enddo

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

end subroutine InitCommonVerifyCoupler

! ************************************************************************** !

subroutine InitCommonReadRegionFiles(realization)
  ! 
  ! Reads in grid cell ids stored in files
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/03/08
  ! 

  use Realization_class
  use Region_module
  use HDF5_module

  implicit none

  type(realization_type) :: realization
  
  type(region_type), pointer :: region
 
  
  region => realization%region_list%first
  do 
    if (.not.associated(region)) exit
    if (len_trim(region%filename) > 1) then
      if (index(region%filename,'.h5') > 0) then
        if (region%grid_type == STRUCTURED_GRID_REGION) then
          call HDF5ReadRegionFromFile(realization,region,region%filename)
        else
          !geh: Do not skip this subouroutine if PETSC_HAVE_HDF5 is not 
          !     defined.  The subroutine prints an error message if not defined
          !     informing the user of the error.  If you skip the subroutine,
          !     no error message is printed and the user is unaware of the
          !     region not being read.
          call HDF5ReadUnstructuredGridRegionFromFile(realization%option,region, &
                                                      region%filename)
        endif
      else if (index(region%filename,'.ss') > 0) then
        region%def_type = DEFINED_BY_SIDESET_UGRID
        region%sideset => RegionCreateSideset()
        call RegionReadFromFile(region%sideset,region%filename, &
                                realization%option)
      else if (index(region%filename,'.ex') > 0) then
        region%def_type = DEFINED_BY_FACE_UGRID_EXP
        call RegionReadFromFile(region%explicit_faceset,region%cell_ids, &
                                region%filename,realization%option)
        region%num_cells = size(region%cell_ids)
      else
        call RegionReadFromFile(region,realization%option, &
                                region%filename)
      endif
    endif
    region => region%next
  enddo

end subroutine InitCommonReadRegionFiles

! ************************************************************************** !

subroutine readVectorFromFile(realization,vector,filename,vector_type)
  ! 
  ! Reads data from a file into an associated vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/08
  ! 

  use Realization_class
  use Discretization_module
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Logging_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  Vec :: vector
  character(len=MAXWORDLENGTH) :: filename
  PetscInt :: vector_type
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr
  PetscInt :: count, read_count, i
  PetscInt :: flag
  PetscInt, pointer :: indices(:)
  PetscReal, pointer :: values(:)
  PetscInt, parameter :: block_size = 10000
  Vec :: natural_vec, global_vec

  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (index(filename,'.h5') > 0) then
    ! to be taken care of later in readPermeabilitiesFromFile()
  else
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      option%io_buffer = 'File: ' // trim(filename) // ' not found.'
      call printErrMsg(option)
    endif
    allocate(values(block_size))
    allocate(indices(block_size))
    call DiscretizationCreateVector(discretization,ONEDOF,natural_vec, &
                                    NATURAL,option)
    count = 0
    do
      if (count >= grid%nmax) exit
      read_count = min(block_size,grid%nmax-count)
      do i=1,read_count
        indices(i) = count+i-1 ! zero-based indexing
      enddo
      ierr = 0
      if (option%myrank == option%io_rank) &
        read(fid,*,iostat=ierr) values(1:read_count)
      flag = ierr
      call MPI_Bcast(flag,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                     option%mycomm,ierr)      
      if (flag /= 0) then
        option%io_buffer = 'Insufficent data in file: ' // filename
        call printErrMsg(option)
      endif
      if (option%myrank == option%io_rank) then
        call VecSetValues(natural_vec,read_count,indices,values,INSERT_VALUES, &
                          ierr);CHKERRQ(ierr)
      endif
      count = count + read_count
    enddo
    call MPI_Bcast(count,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                   option%mycomm,ierr)      
    if (count /= grid%nmax) then
      write(option%io_buffer,'("Number of data in file (",i8, &
      & ") does not match size of vector (",i8,")")') count, grid%nlmax
      call printErrMsg(option)
    endif
    close(fid)
    deallocate(values)
    nullify(values)
    deallocate(indices)
    nullify(indices)
    call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
    select case(vector_type)
      case(LOCAL)
        call DiscretizationCreateVector(discretization,ONEDOF,global_vec, &
                                        GLOBAL,option)        
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                           global_vec,ONEDOF)  
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         vector,ONEDOF)
        call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
      case(GLOBAL)
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                           vector,ONEDOF) 
    end select 
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine readVectorFromFile

! ************************************************************************** !

subroutine InitCommonCreateIOGroups(option)
  ! 
  ! Create sub-communicators that are used in initialization
  ! and output HDF5 routines.
  ! 
  ! Author: Vamsi Sripathi
  ! Date: 07/14/09
  ! 

  use Option_module
  use Logging_module

#if defined(SCORPIO)
  use hdf5
#endif

  implicit none

  type(option_type) :: option
  PetscErrorCode :: ierr

#if defined(SCORPIO)

  PetscMPIInt :: numiogroups

  call PetscLogEventBegin(logging%event_create_iogroups,ierr);CHKERRQ(ierr)

  ! Initialize HDF interface to define global constants  
  call h5open_f(ierr)

  if (option%hdf5_read_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_READ_GROUP_SIZE & 
            & in the input file (pflotran.in) is either not set or &
            & its value is less than or equal to ZERO. &
            & HDF5_READ_GROUP_SIZE =  ",i6)') &
             option%hdf5_read_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let one process read and broadcast to everyone
    option%hdf5_read_group_size = option%mycommsize
  endif         
 
  if (option%hdf5_write_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_WRITE_GROUP_SIZE & 
            &in the input file (pflotran.in) is either not set or &
            &its value is less than or equal to ZERO. &
            &HDF5_WRITE_GROUP_SIZE =  ",i6)') &
             option%hdf5_write_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let everyone write separately 
    option%hdf5_write_group_size = 1
  endif                    

  ! create read IO groups
  numiogroups = option%mycommsize/option%hdf5_read_group_size
  call scorpio_iogroup_init(numiogroups, option%mycomm, option%ioread_group_id, ierr)

  if ( option%hdf5_read_group_size == option%hdf5_write_group_size ) then
    ! reuse read_group to use for writing too as both groups are same size
    option%iowrite_group_id = option%ioread_group_id
  else   
      ! create write IO groups
      numiogroups = option%mycommsize/option%hdf5_write_group_size
      call scorpio_iogroup_init(numiogroups, option%mycomm, option%iowrite_group_id, ierr)
  end if

    write(option%io_buffer, '(" Read group id :  ", i6)') option%ioread_group_id
    call printMsg(option)      
    write(option%io_buffer, '(" Write group id :  ", i6)') option%iowrite_group_id
    call printMsg(option)      
  call PetscLogEventEnd(logging%event_create_iogroups,ierr);CHKERRQ(ierr)
#endif
! SCORPIO
 
end subroutine InitCommonCreateIOGroups

! ************************************************************************** !

subroutine InitCommonPrintPFLOTRANHeader(option,fid)
  ! 
  ! Initializes pflotran
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Option_module
  
  implicit none
  
  PetscInt :: fid
  
  type(option_type) :: option
  
  write(fid,'(" PFLOTRAN Header")') 
  
end subroutine InitCommonPrintPFLOTRANHeader

! ************************************************************************** !

subroutine InitCommonReadVelocityField(realization)
  ! 
  ! Reads fluxes in for transport with no flow.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/13
  ! 

  use Realization_class
  use Patch_module
  use Field_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Discretization_module
  use HDF5_module

  implicit none
  
  type(realization_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: idir, iconn, sum_connection
  PetscInt :: ghosted_id_up, local_id
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_loc_p(:)
  PetscReal, pointer :: vec_p(:)
  type(coupler_type), pointer :: boundary_condition  
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization
  
  filename = realization%nonuniform_velocity_filename

  group_name = ''
  do idir = 1, 3
    select case(idir)
      case(1)
        dataset_name = 'Internal Velocity X'
      case(2)
        dataset_name = 'Internal Velocity Y'
      case(3)
        dataset_name = 'Internal Velocity Z'
    end select
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call DiscretizationGlobalToLocal(discretization,field%work,field%work_loc, &
                                     ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_loc_p,ierr);CHKERRQ(ierr)
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0  
    do 
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        ghosted_id_up = cur_connection_set%id_up(iconn)
        if (cur_connection_set%dist(idir,iconn) > 0.9d0) then
          patch%internal_velocities(1,sum_connection) = vec_loc_p(ghosted_id_up)
        endif
      enddo
      cur_connection_set => cur_connection_set%next
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  enddo
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    dataset_name = boundary_condition%name
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      patch%boundary_velocities(1,sum_connection) = vec_p(local_id)
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine InitCommonReadVelocityField

end module Init_Common_module
