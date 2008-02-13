module Init_module

  implicit none

  private

#include "definitions.h"

#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscsnes.h"


  public :: PflowInit

contains

! ************************************************************************** !
!
! PflowInit: Initializes a pflow grid object
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine PflowInit(simulation,filename)

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Realization_module
  use Material_module
  use Timestepper_module
  use Field_module
  use Connection_module
  use Coupler_module
  use General_Grid_module
  use Debug_module
  use Convergence_module
  
  use MPHASE_module
  use Richards_Lite_module
  use Richards_Analytical_module
  use THC_module

  use Convergence_module
  use Utility_module
    
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  type(stepper_type), pointer :: stepper
  type(solver_type), pointer :: solver
  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(pflow_debug_type), pointer :: debug
  PetscInt :: mcomp, mphas
  PetscInt :: temp_int
  PetscTruth :: iflag

  PetscReal :: alpha, maxstep, steptol

  PetscReal, pointer :: phis_p(:)
                       
  PetscErrorCode :: ierr
  
  ! set pointers to objects
  stepper => simulation%stepper
  solver => stepper%solver
  realization => simulation%realization
  option => realization%option
  field => realization%field
  debug => realization%debug
  
  ! read MODE,GRID,PROC,COMP,PHAS cards
  call readRequiredCardsFromInput(realization,filename,mcomp,mphas)
  grid => realization%grid

  ! set the operational mode (e.g. RICHARDS_MODE, MPH_MODE, etc)
  call setMode(option,mcomp,mphas)
  ! process command line options
  call OptionCheckCommandLine(option)

! check number of dofs and phases
  iflag = PETSC_FALSE
  select case(option%imode)
    case(THC_MODE)
      if (option%ndof .ne. 3 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
    case(MPH_MODE,RICHARDS_MODE)
      if (option%ndof .ne. (option%nspec+1)) iflag = PETSC_TRUE
    case(RICHARDS_LITE_MODE)
      if (option%ndof /= 1 .and. option%nphase /= 1 .and. option%nspec /= 1) &
        iflag = PETSC_TRUE
    case default
      if (option%ndof .ne. 1 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
  end select
  
  if (iflag == PETSC_TRUE) then
    write(*,*) 'Specified number of dofs or phases not correct-stop: ', &
               trim(option%mode), 'ndof= ',option%ndof,' nph= ', &
               option%nphase
    stop
  endif

  call GridCreateDMs(grid,option)
  
 !-----------------------------------------------------------------------
 ! Create the vectors with parallel layout corresponding to the DM's,
 ! and, for vectors that need to be ghosted, create the corresponding
 ! ghosted vectors.
 !-----------------------------------------------------------------------

  ! 1 degree of freedom, global
  call GridCreateVector(grid,ONEDOF,field%porosity0,GLOBAL)
  call VecDuplicate(field%porosity0, grid%volume, ierr)
  call VecDuplicate(field%porosity0, field%perm0_xx, ierr)
  call VecDuplicate(field%porosity0, field%perm0_yy, ierr)
  call VecDuplicate(field%porosity0, field%perm0_zz, ierr)
  call VecDuplicate(field%porosity0, field%perm_pow, ierr)
      
  ! 1 degree of freedom, local
  call GridCreateVector(grid,ONEDOF,field%porosity_loc,LOCAL)
  call VecDuplicate(field%porosity_loc, field%tor_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%ithrm_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%icap_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%iphas_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%iphas_old_loc, ierr)
  
  call VecDuplicate(field%porosity_loc, field%perm_xx_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%perm_yy_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%perm_zz_loc, ierr)

  if (associated(grid%structured_grid)) then
    call VecDuplicate(field%porosity0, grid%structured_grid%dx, ierr)
    call VecDuplicate(field%porosity0, grid%structured_grid%dy, ierr)
    call VecDuplicate(field%porosity0, grid%structured_grid%dz, ierr)

    call VecDuplicate(field%porosity_loc, grid%structured_grid%dx_loc, ierr)
    call VecDuplicate(field%porosity_loc, grid%structured_grid%dy_loc, ierr)
    call VecDuplicate(field%porosity_loc, grid%structured_grid%dz_loc, ierr)
  endif

  ! ndof degrees of freedom, global
  call GridCreateVector(grid,NDOF, field%xx, GLOBAL)
  call VecDuplicate(field%xx, field%yy, ierr)
  call VecDuplicate(field%xx, field%dxx, ierr)
  call VecDuplicate(field%xx, field%r, ierr)
  call VecDuplicate(field%xx, field%accum, ierr)

  ! ndof degrees of freedom, local
  call GridCreateVector(grid,NDOF, field%xx_loc, LOCAL)
     
  ! set up nG2L, NL2G, etc.
  call GridMapIndices(grid)
  
  ! read in the remainder of the input file
  call readInput(simulation,filename)

  if (option%iblkfmt == 0) then
    select case(option%imode)
      case(MPH_MODE,RICHARDS_MODE)
        call printErrMsg(option,&
                         'AIJ matrix not supported for current mode: '// &
                         option%mode)
    end select
  endif

  call GridComputeSpacing(grid)
  call GridComputeCoordinates(grid,option)
  call GridComputeVolumes(grid,option)

  ! read any regions provided in external files
  call readRegionFiles(realization)
  ! clip regions and set up boundary connectivity, distance  
  call GridLocalizeRegions(realization%regions,realization%grid,realization%option)

  ! set up internal connectivity, distance, etc.
  call GridComputeInternalConnect(grid,option)
  call RealizationProcessCouplers(realization)

  ! connectivity between initial conditions, boundary conditions, srcs/sinks, etc and grid
  call GridComputeCouplerConnections(grid,option,realization%initial_conditions%first)
  call GridComputeCouplerConnections(grid,option,realization%boundary_conditions%first)
  call GridComputeCouplerConnections(grid,option,realization%source_sinks%first)
                                
  call assignMaterialPropToRegions(realization)
  call RealizationInitCouplerAuxVars(realization,realization%initial_conditions)
  call RealizationInitCouplerAuxVars(realization,realization%boundary_conditions)
  call assignInitialConditions(realization)

  ! should we still support this
  if (option%use_generalized_grid) then 
    if (option%myrank == 0) print *, 'Reading structured grid from hdf5' 
    if (.not.associated(field%imat)) &
      allocate(field%imat(grid%ngmax))  ! allocate material id array
    call ReadStructuredGridHDF5(realization)
  endif
  
  allocate(realization%field%internal_velocities(option%nphase, &
           ConnectionGetNumberInList(realization%grid%internal_connection_list)))
  realization%field%internal_velocities = 0.d0
  temp_int = CouplerGetNumConnectionsInList(realization%boundary_conditions)
  allocate(realization%field%boundary_velocities(option%nphase,temp_int)) 
  realization%field%boundary_velocities = 0.d0          

  if (option%myrank == 0) then
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    if (grid%igrid == STRUCTURED) then
      write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
        option%commsize,grid%structured_grid%npx,grid%structured_grid%npy, &
        grid%structured_grid%npz
    endif
    write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
      option%ndof,option%nphase
    select case(option%imode)
      case(THC_MODE)
        write(*,'(" mode = THC: p, T, C")')
      case(MPH_MODE)
        write(*,'(" mode = MPH: p, T, s/C")')
      case(RICHARDS_MODE)
        write(*,'(" mode = Richards: p, T, s/C")')
      case(RICHARDS_LITE_MODE)
        write(*,'(" mode = Richards: p")')      
    end select
  endif
  
  call SolverCreateSNES(solver)  

  if (option%use_matrix_free == PETSC_TRUE) then
  
    option%ideriv = 0
  
    if (option%myrank == 0) write(*,'(" Using matrix-free Newton-Krylov")')

    select case(option%imode)
      case(THC_MODE,MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
        call MatCreateMFFD(solver%snes,field%xx,solver%J,ierr)
    end select
        
  else

    option%ideriv = 1
  
    call GridCreateJacobian(grid,solver%J,option)
  
  endif

  if (option%myrank == 0) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')

  select case(option%imode)
    case(THC_MODE)
      call SNESSetFunction(solver%snes,field%r,THCResidual,realization,ierr)
      call SNESSetJacobian(solver%snes, solver%J, solver%J, THCJacobian, &
                           realization, ierr)
    case(RICHARDS_MODE)
      call SNESSetFunction(solver%snes,field%r,RichardsAnalyticalResidual,realization,ierr)
      call SNESSetJacobian(solver%snes, solver%J, solver%J, RichardsAnalyticalJacobian, &
                           realization, ierr)
    case(RICHARDS_LITE_MODE)
      call SNESSetFunction(solver%snes,field%r,RichardsLiteResidual,realization,ierr)
      call SNESSetJacobian(solver%snes, solver%J, solver%J, RichardsLiteJacobian, &
                           realization, ierr)
    case(MPH_MODE)
      call SNESSetFunction(solver%snes,field%r,MPHASEResidual,realization,ierr)
      call SNESSetJacobian(solver%snes, solver%J, solver%J, MPHASEJacobian, &
                           realization, ierr)
  end select

  call SolverSetSNESOptions(solver,option)

  ! shell for custom convergence test.  The default SNES convergence test  
  ! is call within this function. 
  stepper%convergence_context => ConvergenceContextCreate(solver,option)
  call SNESSetConvergenceTest(solver%snes,ConvergenceTest, &
                              stepper%convergence_context,ierr) 

  if (option%myrank == 0) write(*,'("  Finished setting up of SNES ")')

  if (option%myrank == 0) write(*,'("  Finished setting up of INIT ")')
         
  ! move each case to its respective module and just call ModeSetup (e.g. RichardsSetup)
  select case(option%imode)
    case(RICHARDS_MODE)
      call RichardsSetup(realization)
    case(RICHARDS_LITE_MODE)
      call RichardsLiteSetup(realization)
    case(MPH_MODE)
      call MphaseSetup(realization)
    case(THC_MODE)
      call THCSetup(realization)
  end select  

  if (option%myrank == 0) write(*,'("  Finished setting up ")')

  if (debug%print_couplers) then
    call verifyCouplers(realization,realization%initial_conditions)
    call verifyCouplers(realization,realization%boundary_conditions)
    call verifyCouplers(realization,realization%source_sinks)
  endif
  
end subroutine PflowInit

! ************************************************************************** !
!
! readRequiredCardsFromInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine readRequiredCardsFromInput(realization,filename,mcomp,mphas)

  use Option_module
  use Grid_module
  use Fileio_module
  use Realization_module
  
  implicit none

  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename
  PetscInt :: mcomp, mphas, idum

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  PetscInt :: igeom
  
  option => realization%option
  
  open(IUNIT1, file=filename, action="read", status="old") 
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  mphas=0
  mcomp = 0

! Read in select required cards
!.........................................................................

  ! MODE information
  string = "MODE"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(option%myrank,string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! read in keyword 
  call fiReadWord(string,option%mode,.true.,ierr)
  call fiErrorMsg(option%myrank,'mode','mode',ierr)

!.........................................................................

  ! GRID information
  string = "GRID"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(option%myrank,string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! key off igeom for structured vs unstructured 
  call fiReadInt(string,igeom,ierr)
  call fiDefaultMsg(option%myrank,'igeom',ierr)
  
  realization%grid => GridCreate(igeom) 
  grid => realization%grid

  if (grid%igrid == STRUCTURED) then ! structured
    call fiReadInt(string,grid%structured_grid%nx,ierr)
    call fiDefaultMsg(option%myrank,'nx',ierr)
    
    call fiReadInt(string,grid%structured_grid%ny,ierr)
    call fiDefaultMsg(option%myrank,'ny',ierr)
    
    call fiReadInt(string,grid%structured_grid%nz,ierr)
    call fiDefaultMsg(option%myrank,'nz',ierr)
    
    grid%structured_grid%nxy = grid%structured_grid%nx*grid%structured_grid%ny
    grid%structured_grid%nmax = grid%structured_grid%nxy * &
                                grid%structured_grid%nz
    grid%nmax = grid%structured_grid%nmax
  else ! unstructured
  endif
      
  call fiReadInt(string,option%nphase,ierr)
  call fiDefaultMsg(option%myrank,'nphase',ierr)

  call fiReadInt(string,option%nspec,ierr)
  call fiDefaultMsg(option%myrank,'nspec',ierr)

  call fiReadInt(string,option%npricomp,ierr)
  call fiDefaultMsg(option%myrank,'npricomp',ierr)

  call fiReadInt(string,option%ndof,ierr)
  call fiDefaultMsg(option%myrank,'ndof',ierr)
      
  call fiReadInt(string,idum,ierr)
  call fiDefaultMsg(option%myrank,'idcdm',ierr)

  call fiReadInt(string,option%itable,ierr)
  call fiDefaultMsg(option%myrank,'itable',ierr)

  if (option%myrank==0) then
    idum = 0
    if (grid%igrid == STRUCTURED) then
      write(IUNIT2,'(/," *GRID ",/, &
        &"  igeom   = ",i4,/, &
        &"  nx      = ",i4,/, &
        &"  ny      = ",i4,/, &
        &"  nz      = ",i4,/, &
        &"  nphase  = ",i4,/, &
        &"  ndof    = ",i4,/, &
        &"  idcdm   = ",i4,/, &
        &"  itable  = ",i4    &
        &   )') grid%igeom, grid%structured_grid%nx, &
                grid%structured_grid%ny, grid%structured_grid%nz, &
                option%nphase, option%ndof, idum, option%itable
    else
    endif
  endif

!.........................................................................

  if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
    
    ! PROC information
    string = "PROC"
    call fiFindStringInFile(IUNIT1,string,ierr)

    if (ierr == 0) then

      ! strip card from front of string
      call fiReadWord(string,word,.false.,ierr)
      call fiReadInt(string,grid%structured_grid%npx,ierr)
      call fiDefaultMsg(option%myrank,'npx',ierr)
      call fiReadInt(string,grid%structured_grid%npy,ierr)
      call fiDefaultMsg(option%myrank,'npy',ierr)
      call fiReadInt(string,grid%structured_grid%npz,ierr)
      call fiDefaultMsg(option%myrank,'npz',ierr)
 
      if (option%myrank == 0) &
        write(IUNIT2,'(/," *PROC",/, &
          & "  npx   = ",3x,i4,/, &
          & "  npy   = ",3x,i4,/, &
          & "  npz   = ",3x,i4)') grid%structured_grid%npx, &
            grid%structured_grid%npy, grid%structured_grid%npz
  
      if (option%commsize /= grid%structured_grid%npx * &
                             grid%structured_grid%npy * &
                             grid%structured_grid%npz) then
        if (option%myrank==0) &
          write(*,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%npx*grid%structured_grid%npy* &
                       grid%structured_grid%npz,' commsize = ',option%commsize
        stop
      endif
    endif
  endif
  
!.........................................................................

  ! COMP information
  string = "COMP"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr == 0) then

    mcomp = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg(option%myrank,'COMP',ierr)
      
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg(option%myrank,'namcx','GAS',ierr)
        
      call fiWordToUpper(name) 
      select case(name(1:len_trim(name)))
        case('H2O')
            mcomp = mcomp +1
        case('TRACER')     
            mcomp = mcomp +4
        case('CO2')
            mcomp = mcomp +2   
        case('C10H22')
            mcomp = mcomp +8
        case('AIR')              
            mcomp = mcomp +16
        case('ENG')
            mcomp = mcomp +32
        case default               
           print *,' Wrong comp name::  Stop:' 
           stop
      end select 
    enddo
  endif          

!....................................................................

  ! COMP information
  string = "PHAS"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr == 0) then

    mphas = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg(option%myrank,'phase',ierr)
    
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg(option%myrank,'namcx','phase',ierr)
        
      call fiWordToUpper(name) 
         
      select case(name(1:len_trim(name)))
        case('ROCK')
          mphas = mphas +1
        case('H2O')     
          mphas = mphas +2
        case('CO2')
          mphas = mphas +4   
        case('GAS')
          mphas = mphas +8
        case('NAPL1')              
          mphas = mphas +16
        case('NAPL2')
          mphas = mphas +32
        case default               
          print *,' Wrong phase name::  Stop:' 
          stop
      end select 
    enddo
  endif
    
end subroutine readRequiredCardsFromInput

! ************************************************************************** !
!
! readInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readInput(simulation,filename)

  use Simulation_module
  use Option_module
  use Field_module
  use Grid_module
  use Structured_Grid_module
  use Solver_module
  use Material_module
  use Fileio_module
  use Realization_module
  use Timestepper_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module
  use Breakthrough_module
  use Waypoint_module
  use Debug_module

  use pckr_module 
  
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
    
  PetscReal, parameter:: fmwnacl = 58.44277D0, fmwh2o  = 18.01534d0
  PetscInt :: i, i1, i2, idum, ireg, isrc, j
  PetscInt :: ibc, ibrk, ir,np
  PetscReal :: rdum

  logical :: continuation_flag
  logical :: periodic_output_flag = .false.
  PetscReal :: periodic_rate = 0.d0
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscInt :: temp_int
  PetscInt :: length 
  PetscInt :: count, id
  
! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
!           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR

  type(region_type), pointer :: region
  type(condition_type), pointer :: condition
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(breakthrough_type), pointer :: breakthrough
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_type), pointer :: material
  type(thermal_property_type), pointer :: thermal_property
  type(saturation_function_type), pointer :: saturation_function

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(solver_type), pointer :: solver
  type(stepper_type), pointer :: stepper
  
  realization => simulation%realization
  grid => realization%grid
  option => realization%option
  field => realization%field
  stepper => simulation%stepper
  solver => stepper%solver

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                              
  rewind(IUNIT1)  
    
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    length = len_trim(word)
    call fiCharsToUpper(word,length)
!    call fiReadCard(word,card,ierr)
    card = trim(word)

    if (option%myrank == 0) print *, 'pflow_read:: ',card

    select case(trim(card))

!....................
      case ('MODE')

!....................
      case ('GRID')
      
!....................
      case ('DEBUG','PFLOW_DEBUG')
        call DebugRead(realization%debug,IUNIT1,option%myrank)
        
!....................
      case ('GENERALIZED_GRID')
        option%use_generalized_grid = .true.
        call fiReadWord(string,option%generalized_grid,.true.,ierr)

!....................
      case ('PROC')
      
!....................
      case ('REGION','REGN')
        region => RegionCreate()
        call fiReadWord(string,region%name,.true.,ierr)
        call fiErrorMsg(option%myrank,'regn','name',ierr) 
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg(option%myrank,'REGN',ierr)
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(option%myrank,'type','REGN', ierr)
        if (fiStringCompare(word,"BLOCK",FIVE_INTEGER)) then ! block region
          if (grid%igrid /= STRUCTURED) then
            call printErrMsg(option,"BLOCK region not supported for &
                             &unstructured grid")
          endif
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'REGN',ierr)
          call fiReadInt(string,region%i1,ierr) 
          call fiErrorMsg(option%myrank,'i1','REGN', ierr)
          call fiReadInt(string,region%i2,ierr)
          call fiErrorMsg(option%myrank,'i2','REGN', ierr)
          call fiReadInt(string,region%j1,ierr)
          call fiErrorMsg(option%myrank,'j1','REGN', ierr)
          call fiReadInt(string,region%j2,ierr)
          call fiErrorMsg(option%myrank,'j2','REGN', ierr)
          call fiReadInt(string,region%k1,ierr)
          call fiErrorMsg(option%myrank,'k1','REGN', ierr)
          call fiReadInt(string,region%k2,ierr)
          call fiErrorMsg(option%myrank,'k2','REGN', ierr)
        else if (fiStringCompare(word,"LIST",FOUR_INTEGER)) then
          word = ""
          call fiReadWord(string,word,.true.,ierr)
          if (fiStringCompare(word,"file",FOUR_INTEGER)) then
!            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg(option%myrank,'REGN',ierr)
            call fiReadWord(string,word,.true.,ierr)
            call fiErrorMsg(option%myrank,'filename','REGN', ierr)
            region%filename = word
            ! this file will be read later in readRegionFiles()          
          else
            call RegionReadFromFile(region,IUNIT1)
          endif            
        else
          call printErrMsg(option,"REGION type not recognized")
        endif
        call RegionAddToList(region,realization%regions)      

!....................
      case ('CONDITION','COND')
        condition => ConditionCreate(option)
        call fiReadWord(string,condition%name,.true.,ierr)
        call fiErrorMsg(option%myrank,'cond','name',ierr) 
        call ConditionRead(condition,option,IUNIT1)
        call ConditionAddToList(condition,realization%conditions)
      
!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1,option)
        call CouplerAddToList(coupler,realization%boundary_conditions)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1,option)
        call CouplerAddToList(coupler,realization%initial_conditions)
      
!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,IUNIT1,option)
        call StrataAddToList(strata,realization%strata)
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1,option)
        call CouplerAddToList(coupler,realization%source_sinks)
      
!.....................
      case ('COMP') 
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo
!.....................
      case ('PHAS')
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'PHASE',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo 
      
!....................

      case ('COUP')

        call printWrnMsg(option,"COUP not currently supported")
        call fiReadStringErrorMsg(option%myrank,'COUP',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'isync',ierr)

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') idum

!....................

      case ('GRAV','GRAVITY')

        call fiReadStringErrorMsg(option%myrank,'GRAV',ierr)

        call fiReadDouble(string,temp_real,ierr)
        if (ierr /= 0) then
          call fiDefaultMsg(option%myrank,'gravity',ierr)
        else
          call fiReadDouble(string,option%gravity(2),ierr)
          if (ierr /= 0) then
            option%gravity(:) = 0.d0
            option%gravity(3) = temp_real
          else
            option%gravity(1) = temp_real
            call fiReadDouble(string,option%gravity(3),ierr)
          endif
        endif

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,3pe12.4 &
            & )') option%gravity(1:3)

!....................

      case ('HDF5')
        realization%output_option%print_hdf5 = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          length = len_trim(word)
          call fiCharsToUpper(word,length)
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_hdf5_velocities = .true.
            case('FLUX')
              realization%output_option%print_hdf5_flux_velocities = .true.
            case default
          end select
            
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *HDF5",10x,i1,/)') realization%output_option%print_hdf5

!.....................
      case ('INVERT_Z','INVERTZ')
        if (associated(grid%structured_grid)) then
          grid%structured_grid%invert_z_axis = .true.
          option%gravity(3) = -1.d0*option%gravity(3)
        endif
      
!....................

      case ('TECP')
        realization%output_option%print_tecplot = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          length = len_trim(word)
          call fiCharsToUpper(word,length)
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_tecplot_velocities = .true.
            case('FLUX')
              realization%output_option%print_tecplot_flux_velocities = .true.
            case default
          end select
          
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *TECP",10x,i1,/)') realization%output_option%print_tecplot

!....................


      case ('OPTS')

        call fiReadStringErrorMsg(option%myrank,'OPTS',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'write_init',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'iprint',ierr)

        call fiReadInt(string,option%imod,ierr)
        call fiDefaultMsg(option%myrank,'mod',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'itecplot',ierr)

        call fiReadInt(string,option%iblkfmt,ierr)
        call fiDefaultMsg(option%myrank,'iblkfmt',ierr)

        call fiReadInt(string,stepper%ndtcmx,ierr)
        call fiDefaultMsg(option%myrank,'ndtcmx',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'iran_por',ierr)
  
        call fiReadDouble(string,rdum,ierr)
        call fiDefaultMsg(option%myrank,'ran_fac',ierr)
    
        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'iread_perm',ierr)
    
        call fiReadInt(string,option%iread_geom,ierr)
        call fiDefaultMsg(option%myrank,'iread_geom',ierr)

        idum = 0
        if (option%myrank == 0) &
          write(IUNIT2,'(/," *OPTS",/, &
            & "  write_init = ",3x,i2,/ &
            & "  imod       = ",3x,i2,/, &
            & "  itecplot   = ",3x,i2,/, &
            & "  iblkfmt    = ",3x,i2,/, &
            & "  ndtcmx     = ",3x,i2,/, &
            & "  iran_por   = ",3x,i2,/, &
            & "  ran_fac    = ",3x,1pe12.4,/, &
            & "  iread_perm = ",3x,i2,/, &
            & "  iread_geom = ",3x,i2 &
            & )') idum,option%imod,idum, &
            option%iblkfmt,stepper%ndtcmx,idum,idum,idum,option%iread_geom

!....................

      case ('TOLR')

        call fiReadStringErrorMsg(option%myrank,'TOLR',ierr)

        call fiReadInt(string,stepper%stepmax,ierr)
        call fiDefaultMsg(option%myrank,'stepmax',ierr)
  
        call fiReadInt(string,stepper%iaccel,ierr)
        call fiDefaultMsg(option%myrank,'iaccel',ierr)

        call fiReadInt(string,stepper%newton_max,ierr)
        call fiDefaultMsg(option%myrank,'newton_max',ierr)

        call fiReadInt(string,stepper%icut_max,ierr)
        call fiDefaultMsg(option%myrank,'icut_max',ierr)

        call fiReadDouble(string,option%dpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dpmxe',ierr)

        call fiReadDouble(string,option%dtmpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dtmpmxe',ierr)
  
        call fiReadDouble(string,option%dcmxe,ierr)
        call fiDefaultMsg(option%myrank,'dcmxe',ierr)

        call fiReadDouble(string,option%dsmxe,ierr)
        call fiDefaultMsg(option%myrank,'dsmxe',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
          &"  flowsteps  = ",i6,/,      &
          &"  iaccel     = ",i3,/,      &
          &"  newtmx     = ",i3,/,      &
          &"  icutmx     = ",i3,/,      &
          &"  dpmxe      = ",1pe12.4,/, &
          &"  dtmpmxe    = ",1pe12.4,/, &
          &     "  dcmxe      = ",1pe12.4,/, &
          &"  dsmxe      = ",1pe12.4)') &
! For commented-out lines to work with the Sun f95 compiler, we have to 
! terminate the string in the line above; otherwise, the compiler tries to
! include the commented-out line as part of the continued string.
          stepper%stepmax,stepper%iaccel,stepper%newton_max,stepper%icut_max, &
          option%dpmxe,option%dtmpmxe,option%dcmxe, option%dsmxe

!....................

      case ('DXYZ')
      
        if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
          call StructuredGridReadDXYZ(grid%structured_grid,option)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "DXYZ" not supported for unstructured grid'
            stop
        endif

!....................

      case('ORIG','ORIGIN')
        call fiReadDouble(string,grid%origin(X_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'X direction','Origin',ierr)
        call fiReadDouble(string,grid%origin(Y_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'Y direction','Origin',ierr)
        call fiReadDouble(string,grid%origin(Z_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'Z direction','Origin',ierr)
        
!....................

      case('RAD0')
    
        if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
          call fiReadDouble(string,grid%structured_grid%Radius_0,ierr)
          call fiDefaultMsg(option%myrank,'R_0',ierr)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "RAD0" not supported for unstructured grid'
            stop
        endif


      case ('DIFF')

        call fiReadStringErrorMsg(option%myrank,'DIFF',ierr)

        call fiReadDouble(string,option%difaq,ierr)
        call fiDefaultMsg(option%myrank,'difaq',ierr)

        call fiReadDouble(string,option%delhaq,ierr)
        call fiDefaultMsg(option%myrank,'delhaq',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          option%difaq,option%delhaq

!....................

      case ('RCTR')

        call printErrMsg(option,"RCTR currently out of date.  Needs to be reimplemented")
#if 0
        call fiReadStringErrorMsg(option%myrank,'RCTR',ierr)

        call fiReadInt(string,option%ityprxn,ierr)
        call fiDefaultMsg(option%myrank,'ityprxn',ierr)

        call fiReadDouble(string,option%rk,ierr)
        call fiDefaultMsg(option%myrank,'rk',ierr)

        call fiReadDouble(string,option%phis0,ierr)
        call fiDefaultMsg(option%myrank,'phis0',ierr)

        call fiReadDouble(string,option%areas0,ierr)
        call fiDefaultMsg(option%myrank,'areas0',ierr)

        call fiReadDouble(string,option%pwrsrf,ierr)
        call fiDefaultMsg(option%myrank,'pwrsrf',ierr)

        call fiReadDouble(string,option%vbars,ierr)
        call fiDefaultMsg(option%myrank,'vbars',ierr)

        call fiReadDouble(string,option%ceq,ierr)
        call fiDefaultMsg(option%myrank,'ceq',ierr)

        call fiReadDouble(string,option%delHs,ierr)
        call fiDefaultMsg(option%myrank,'delHs',ierr)

        call fiReadDouble(string,option%delEs,ierr)
        call fiDefaultMsg(option%myrank,'delEs',ierr)

        call fiReadDouble(string,option%wfmts,ierr)
        call fiDefaultMsg(option%myrank,'wfmts',ierr)

        if (option%myrank == 0) &
        write(IUNIT2,'(/," *RCTR",/, &
          & "  ityp   = ",3x,i3,/, &
          & "  rk     = ",3x,1pe12.4," [mol/cm^2/s]",/, &
          & "  phis0  = ",3x,1pe12.4," [-]",/, &
          & "  areas0 = ",3x,1pe12.4," [1/cm]",/, &
          & "  pwrsrf = ",3x,1pe12.4," [-]",/, &
          & "  vbars  = ",3x,1pe12.4," [cm^3/mol]",/, &
          & "  ceq    = ",3x,1pe12.4," [mol/L]",/, &
          & "  delHs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  delEs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  wfmts  = ",3x,1pe12.4," [g/mol]" &
          & )') option%ityprxn,option%rk,option%phis0,option%areas0,option%pwrsrf, &
          option%vbars,option%ceq,option%delHs,option%delEs,option%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
        option%rk = option%rk * option%areas0 * 1.d3
        option%vbars = option%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
        option%delHs = option%delHs * option%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!        option%delHs = option%delHs * option%scale ! convert J/kmol -> MJ/kmol
#endif
!....................

      case ('RADN')

        call fiReadStringErrorMsg(option%myrank,'RADN',ierr)

        call fiReadDouble(string,option%ret,ierr)
        call fiDefaultMsg(option%myrank,'ret',ierr)

        call fiReadDouble(string,option%fc,ierr)
        call fiDefaultMsg(option%myrank,'fc',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          option%ret,option%fc

!....................


      case ('PHAR')
        call printWrnMsg(option,"PHAR currently out of date.  Needs to be reimplemented")
#if 0
        call fiReadStringErrorMsg(option%myrank,'PHAR',ierr)

        call fiReadDouble(string,option%qu_kin,ierr)
        call fiDefaultMsg(option%myrank,'TransReaction',ierr)
        if (option%myrank==0) write(IUNIT2,'(/," *PHAR ",1pe12.4)')option%qu_kin
        option%yh2o_in_co2 = 0.d0
        if (option%qu_kin > 0.d0) option%yh2o_in_co2 = 1.d-2 ! check this number!
#endif     
!......................

      case('RICH')
        call fiReadStringErrorMsg(option%myrank,'RICH',ierr)
        call fiReadDouble(string,option%pref,ierr)
        call fiDefaultMsg(option%myrank,'Ref. Pressure',ierr) 

!......................

      case('BRIN','BRINE')
        call fiReadStringErrorMsg(option%myrank,'BRIN',ierr)
        call fiReadDouble(string,option%m_nacl,ierr)
        call fiDefaultMsg(option%myrank,'NaCl Concentration',ierr) 

        call fiReadWord(string,word,.false.,ierr)
        call fiWordToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            option%m_nacl = option%m_nacl /fmwnacl/(1.D0-option%m_nacl)
          case('MOLE')    
            option%m_nacl = option%m_nacl /fmwh2o/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (option%myrank == 0) print *, option%m_nacl
!......................

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call fiReadWord(string,option%restart_file,.true.,ierr)
        call fiErrorMsg(option%myrank,'RESTART','Restart file name',ierr) 

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call fiReadInt(string,option%checkpoint_frequency,ierr)
        call fiErrorMsg(option%myrank,'CHECKPOINT','Checkpoint frequency',ierr) 

!......................

      case ('NUMERICAL_JACOBIAN')
        option%numerical_derivatives = .true.

      case ('INEXACT_NEWTON')
        solver%inexact_newton = .true.

!......................

      case ('NO_PRINT_CONVERGENCE')
        option%print_convergence = PETSC_FALSE

      case ('NO_INF_NORM','NO_INFINITY_NORM')
        option%check_infinity_norm = PETSC_FALSE

      case ('NO_FORCE_ITERATION')
        option%force_at_least_1_iteration = PETSC_FALSE

      case ('PRINT_DETAILED_CONVERGENCE')
        option%print_detailed_convergence = PETSC_TRUE

!....................

      case ('SOLVER')
        call SolverRead(solver,IUNIT1,option%myrank)

!....................

      case ('SOLV')
    
        call fiReadStringErrorMsg(option%myrank,'SOLV',ierr)

!       call fiReadDouble(string,eps,ierr)
!       call fiDefaultMsg(option%myrank,'eps',ierr)

        call fiReadDouble(string,solver%atol,ierr)
        call fiDefaultMsg(option%myrank,'atol_petsc',ierr)

        call fiReadDouble(string,solver%rtol,ierr)
        call fiDefaultMsg(option%myrank,'rtol_petsc',ierr)

        call fiReadDouble(string,solver%stol,ierr)
        call fiDefaultMsg(option%myrank,'stol_petsc',ierr)
      
        solver%dtol=1.D5
!       if (option%use_ksp == 1) then
        call fiReadDouble(string,solver%dtol,ierr)
        call fiDefaultMsg(option%myrank,'dtol_petsc',ierr)
!       endif
   
        call fiReadInt(string,solver%maxit,ierr)
        call fiDefaultMsg(option%myrank,'maxit',ierr)
      
        call fiReadInt(string,solver%maxf,ierr)
        call fiDefaultMsg(option%myrank,'maxf',ierr)
       
        call fiReadInt(string,option%idt_switch,ierr)
        call fiDefaultMsg(option%myrank,'idt',ierr)
        
        solver%inf_tol = solver%atol
        call fiReadDouble(string,solver%inf_tol,ierr)
        call fiDefaultMsg(option%myrank,'inf_tol_pflow',ierr)
 
        if (option%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
          &"  atol_petsc   = ",1pe12.4,/, &
          &"  rtol_petsc   = ",1pe12.4,/, &
          &"  stol_petsc   = ",1pe12.4,/, &
          &"  dtol_petsc   = ",1pe12.4,/, &
          &"  maxit        = ",8x,i5,/, &
          &"  maxf        = ",8x,i5,/, &
          &"  idt         = ",8x,i5 &
          &    )') &
           solver%atol,solver%rtol,solver%stol,solver%dtol,solver%maxit, &
           solver%maxf,option%idt_switch

! The line below is a commented-out portion of the format string above.
! We have to put it here because of the stupid Sun compiler.
!    &"  eps          = ",1pe12.4,/, &

!....................

      case ('THRM','THERMAL_PROPERTY','THERMAL_PROPERTIES')

        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'THRM',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
       
          count = count + 1
          thermal_property => ThermalPropertyCreate()
      
          call fiReadInt(string,thermal_property%id,ierr)
          call fiErrorMsg(option%myrank,'id','THRM', ierr)

          call fiReadDouble(string,thermal_property%rock_density,ierr)
          call fiErrorMsg(option%myrank,'rock density','THRM', ierr)

          call fiReadDouble(string,thermal_property%spec_heat,ierr)
          call fiErrorMsg(option%myrank,'cpr','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_dry,ierr)
          call fiErrorMsg(option%myrank,'ckdry','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_wet,ierr)
          call fiErrorMsg(option%myrank,'ckwet','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%tort_bin_diff,ierr)
          call fiErrorMsg(option%myrank,'tau','THRM', ierr)

          call fiReadDouble(string,thermal_property%vap_air_diff_coef,ierr)
          call fiErrorMsg(option%myrank,'cdiff','THRM', ierr)

          call fiReadDouble(string,thermal_property%exp_binary_diff,ierr)
          call fiErrorMsg(option%myrank,'cexp','THRM', ierr)

        !scale thermal properties
          thermal_property%spec_heat = option%scale * &
                                       thermal_property%spec_heat
          thermal_property%therm_cond_dry = option%scale * &
                                            thermal_property%therm_cond_dry
          thermal_property%therm_cond_wet = option%scale * &
                                            thermal_property%therm_cond_wet
          
          call ThermalAddPropertyToList(thermal_property, &
                                        realization%thermal_properties)
        enddo
        
        ! allocate dynamic arrays holding saturation function information
        allocate(option%rock_density(count))
        allocate(option%cpr(count))
        allocate(option%dencpr(count))
        allocate(option%ckdry(count))
        allocate(option%ckwet(count))
        allocate(option%tau(count))
        allocate(option%cdiff(count))
        allocate(option%cexp(count))
        
        ! fill arrays with values from linked list
        thermal_property => realization%thermal_properties
        do
        
          if (.not.associated(thermal_property)) exit
          
          id = thermal_property%id
          
          if (id > count) then
            call printErrMsg(option,'Thermal property id greater than &
                                    &number of thermal properties')
          endif
                    
          option%rock_density(id) = thermal_property%rock_density
          option%cpr(id) = thermal_property%spec_heat
          option%dencpr(id) = thermal_property%rock_density * &
                              thermal_property%spec_heat
          option%ckdry(id) = thermal_property%therm_cond_dry
          option%ckwet(id) = thermal_property%therm_cond_wet
          option%tau(id) = thermal_property%tort_bin_diff
          option%cdiff(id) = thermal_property%vap_air_diff_coef
          option%cexp(id) = thermal_property%exp_binary_diff
          
          thermal_property => thermal_property%next
          
        enddo
        
        do i=1,count
          if (option%rock_density(i) < 1.d-40) then
            call printErrMsg(option,'Thermal property ids must be numbered &
                             &consecutively from 1 to N')
          endif
        enddo
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *THRM: ",i3)') count
          write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, count
            write(IUNIT2,'(i4,1p7e11.4)') i,option%rock_density(i), &
            option%cpr(i),option%ckdry(i),option%ckwet(i), &
            option%tau(i),option%cdiff(i),option%cexp(i)
          enddo
        endif

!....................

      case ('PCKR','SATURATION_FUNCTION','SATURATION_FUNCTIONS')
      
        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'PCKR',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
       
          count = count + 1
          saturation_function => SaturationFunctionCreate(option)
          
          call fiReadInt(string,saturation_function%id,ierr)
          call fiErrorMsg(option%myrank,'id','PCKR', ierr)
          
          call fiReadInt(string,saturation_function%saturation_function_itype,ierr)
          call fiErrorMsg(option%myrank,'icaptype','PCKR', ierr)
      
          select case(option%imode)
            case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
              do np=1, option%nphase
                call fiReadDouble(string,saturation_function%Sr(np),ierr)
                call fiErrorMsg(option%myrank,'Sr','PCKR', ierr)
              enddo 
            case default
              call fiReadDouble(string,saturation_function%Sr(1),ierr)
              call fiErrorMsg(option%myrank,'Sr','PCKR', ierr)
          end select
        
          call fiReadDouble(string,saturation_function%m,ierr)
          call fiErrorMsg(option%myrank,'pckrm','PCKR', ierr)
          saturation_function%lambda = saturation_function%m

          call fiReadDouble(string,saturation_function%alpha,ierr)
          call fiErrorMsg(option%myrank,'alpha','PCKR', ierr)

          call fiReadDouble(string,saturation_function%pcwmax,ierr)
          call fiErrorMsg(option%myrank,'pcwmax','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%betac,ierr)
          call fiErrorMsg(option%myrank,'pbetac','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%power,ierr)
          call fiErrorMsg(option%myrank,'pwrprm','PCKR', ierr)
          
          call SaturationFunctionAddToList(saturation_function, &
                                           realization%saturation_functions)

        enddo
        
        ! allocate dynamic arrays holding saturation function information
        allocate(option%icaptype(count))
        option%icaptype = 0
  
        select case(option%imode)
          case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
            allocate(option%sir(1:option%nphase,count))
          case default
            allocate(option%swir(count))
        end select
  
        allocate(option%lambda(count))
        allocate(option%alpha(count))
        allocate(option%pckrm(count))
        allocate(option%pcwmax(count))
        allocate(option%pcbetac(count))
        allocate(option%pwrprm(count))

        ! fill arrays with values from linked list
        saturation_function => realization%saturation_functions
        do 
        
          if (.not.associated(saturation_function)) exit
          
          id = saturation_function%id
          
          if (id > count) then
            call printErrMsg(option,'Saturation function id greater than &
                                    &number of saturation functions')
          endif
          
          option%icaptype(id) = saturation_function%saturation_function_itype
          select case(option%imode)
            case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE)
              do i=1,option%nphase
                option%sir(i,id) = saturation_function%Sr(i)
              enddo
            case default
              option%swir(id) = saturation_function%Sr(1)
          end select
          option%lambda(id) = saturation_function%lambda
          option%alpha(id) = saturation_function%alpha
          option%pckrm(id) = saturation_function%m
          option%pcwmax(id) = saturation_function%pcwmax
          option%pcbetac(id) = saturation_function%betac
          option%pwrprm(id) = saturation_function%power
          
          saturation_function => saturation_function%next
          
        enddo
        
        ! check to ensure that all saturation functions were set based on id
        do id = 1,count
          if (option%icaptype(id) == 0) then
            call printErrMsg(option,'Saturation function ids must be numbered &
                               &consecutively from 1 to N')
          endif
        enddo

        if (option%imode == MPH_MODE .or. &
            option%imode == RICHARDS_MODE .or. &
            option%imode == RICHARDS_LITE_MODE) then
          call pckr_init(option%nphase,count,grid%nlmax, &
                         option%icaptype,option%sir, option%pckrm, &
                         option%lambda,option%alpha,option%pcwmax, &
                         option%pcbetac,option%pwrprm)
        endif 

      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *PCKR: ",i3)') ireg
          write(IUNIT2,'("  icp swir    lambda         alpha")')
          do j = 1, count
            i=option%icaptype(j)
            if (option%imode == MPH_MODE .or. &
                option%imode == RICHARDS_MODE .or. &
                option%imode == RICHARDS_LITE_MODE) then
              write(IUNIT2,'(i4,1p8e12.4)') i,(option%sir(np,i),np=1, &
                option%nphase),option%lambda(i),option%alpha(i), &
                option%pcwmax(i),option%pcbetac(i),option%pwrprm(i)
            else
              write(IUNIT2,'(i4,1p7e12.4)') i,option%swir(i), &
                option%lambda(i),option%alpha(i),option%pcwmax(i), &
                option%pcbetac(i),option%pwrprm(i)
            endif
          enddo
        end if

        if (option%imode == MPH_MODE .or. &
            option%imode == RICHARDS_MODE .or. &
            option%imode == RICHARDS_LITE_MODE) then
          deallocate(option%icaptype, option%pckrm, option%lambda, &
                     option%alpha,option%pcwmax, option%pcbetac, &
                     option%pwrprm)
        endif 
 
        call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                           realization%saturation_function_array)
        
!....................
      
      case ('PHIK','MATERIAL','MATERIALS')

        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'PHIK',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit
       
          count = count + 1
          material => MaterialCreate()

          call fiReadWord(string,material%name,.true.,ierr)
          call fiErrorMsg(option%myrank,'name','PHIK', ierr)
                
          call fiReadInt(string,material%id,ierr)
          call fiErrorMsg(option%myrank,'id','PHIK', ierr)
                
          call fiReadInt(string,material%icap,ierr)
          call fiErrorMsg(option%myrank,'icap','PHIK', ierr)
  
          call fiReadInt(string,material%ithrm,ierr)
          call fiErrorMsg(option%myrank,'ithrm','PHIK', ierr)
  
          call fiReadDouble(string,material%porosity,ierr)
          call fiErrorMsg(option%myrank,'por','PHIK', ierr)
          
          call fiReadDouble(string,material%tortuosity,ierr)
          call fiErrorMsg(option%myrank,'tor','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(1,1),ierr)
          call fiErrorMsg(option%myrank,'permx','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(2,2),ierr)
          call fiErrorMsg(option%myrank,'permy','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(3,3),ierr)
          call fiErrorMsg(option%myrank,'permz','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability_pwr,ierr)
          call fiErrorMsg(option%myrank,'permpwr','PHIK', ierr)
          
          material%permeability(1:3,1:3) = material%permeability(1:3,1:3)
          
          call MaterialAddToList(material,realization%materials)
          
        enddo          

        call MaterialConvertListToArray(realization%materials, &
                                        realization%material_array)
                                        
!....................
      
      case ('TIME')

        call fiReadStringErrorMsg(option%myrank,'TIME',ierr)
      
        call fiReadWord(string,word,.false.,ierr)
      
        realization%output_option%tunit = trim(word)

        if (realization%output_option%tunit == 's') then
          realization%output_option%tconv = 1.d0
        else if (realization%output_option%tunit == 'm') then
          realization%output_option%tconv = 60.d0
        else if (realization%output_option%tunit == 'h') then
          realization%output_option%tconv = 60.d0 * 60.d0
        else if (realization%output_option%tunit == 'd') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0
        else if (realization%output_option%tunit == 'mo') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
        else if (realization%output_option%tunit == 'y') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
        else
          if (option%myrank == 0) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') realization%output_option%tunit
          endif
          stop
        endif


        call fiReadWord(string,word,.false.,ierr)
        if (ierr == 0) then
          call fiWordToUpper(word)
          if (fiStringCompare(word,'EVERY',FIVE_INTEGER)) then
            periodic_output_flag = .true.
            call fiReadDouble(string,periodic_rate,ierr)
          endif
        endif

        continuation_flag = .true.
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = .false.
          if (index(string,backslash) > 0) continuation_flag = .true.
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              waypoint%print_output = .true.              
              call WaypointInsertInList(waypoint,stepper%waypoints)
            endif
          enddo
        enddo
        
        ! make last waypoint final
        waypoint%final = .true.
        
        if (periodic_output_flag) then
          temp_real2 = waypoint%time   ! final simulation time
          temp_real = periodic_rate
          do
            if (temp_real > temp_real2) exit
            
            waypoint => WaypointCreate()
            waypoint%time = temp_real
            waypoint%print_output = .true.              
            call WaypointInsertInList(waypoint,stepper%waypoints)
            
            temp_real = temp_real + periodic_rate
          enddo
        endif

!....................

      case ('DTST')

        call fiReadStringErrorMsg(option%myrank,'DTST',ierr)

        call fiReadDouble(string,stepper%dt_min,ierr)
        call fiDefaultMsg(option%myrank,'dt_min',ierr)
            
        continuation_flag = .true.
        temp_int = 0       
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = .false.
          if (index(string,backslash) > 0) continuation_flag = .true.
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              call fiReadDouble(string,waypoint%dt_max,ierr)
              call fiErrorMsg(option%myrank,'dt_max','dtst',ierr)
              if (temp_int == 0) stepper%dt_max = waypoint%dt_max
              call WaypointInsertInList(waypoint,stepper%waypoints)
              temp_int = temp_int + 1
            endif
          enddo
        enddo
        
        option%dt = stepper%dt_min
      
        option%dt = realization%output_option%tconv * option%dt
        stepper%dt_min = realization%output_option%tconv * stepper%dt_min
        stepper%dt_max = realization%output_option%tconv * stepper%dt_max

!....................
      case ('BRK','BREAKTHROUGH')
        breakthrough => BreakthroughCreate()
        call BreakthroughRead(breakthrough,IUNIT1,option)
        call BreakthroughAddToList(breakthrough,realization%breakthrough)
      
!....................
      case('SDST')
        print *, 'SDST needs to be implemented'
        stop
#if 0
! Needs implementation         
        allocate(stepper%steady_eps(option%ndof))
        do j=1,option%ndof
          call fiReadDouble(string,stepper%steady_eps(j),ierr)
          call fiDefaultMsg(option%myrank,'steady tol',ierr)
        enddo
        if (option%myrank==0) write(IUNIT2,'(/," *SDST ",/, &
          &"  dpdt        = ",1pe12.4,/, &
          &"  dtmpdt        = ",1pe12.4,/, &
          &"  dcdt        = ",1pe12.4)') &
          stepper%steady_eps
#endif

!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call fiReadDouble(string,option%wallclock_stop_time,ierr)
        call fiErrorMsg(option%myrank,'stop time','WALLCLOCK_STOP', ierr) 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*3600.d0
      
!....................
      case default
    
        if (option%myrank == 0) then
          print *, "Error reading input file: keyword (", trim(word), &
                   ") not found. Terminating."
        endif
        call PetscFinalize(ierr)
        stop

    end select

  enddo

  close(IUNIT1)
  
end subroutine readInput

! ************************************************************************** !
!
! setMode: Sets the flow mode (richards, vadose, mph, etc.)
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine setMode(option,mcomp,mphas)

  use Option_module
  use Fileio_module

  implicit none 

  type(option_type) :: option
  PetscInt :: mcomp, mphas
  PetscInt :: length
  
  length = len_trim(option%mode)
  call fiCharsToLower(option%mode,length)
  length = len_trim(option%mode)
  if (fiStringCompare(option%mode,"richards",length)) then
    option%imode = RICHARDS_MODE
  else if (fiStringCompare(option%mode,"richards_lite",length)) then
    option%imode = RICHARDS_LITE_MODE
    option%nphase = 1
    option%nspec = 1
    option%ndof = 1
#if 0  
  ! needs to be implemented
  else if (fiStringCompare(option%mode,"MPH",len_trim(option%mode))) then
  else if (fiStringCompare(option%mode,"",#)) then
#endif  
  endif 
  
  if (option%imode /= THC_MODE .and. &
      option%imode /= MPH_MODE .and. &
      option%imode /= RICHARDS_MODE .and. &
      option%imode /= RICHARDS_LITE_MODE) then 

    if (mcomp >0 .and. mphas>0)then
      if (option%imode /= THC_MODE .and. mcomp == 37)then
        option%imode = THC_MODE
        option%mode = 'thc'
        option%nphase = 1; option%ndof =3
      endif
      if (option%imode /= MPH_MODE .and. mcomp == 35)then
        option%imode = MPH_MODE
        option%mode = 'mph'
        option%nphase = 2; option%ndof =3; option%nspec =2 
      endif
      if (option%imode /= RICHARDS_MODE .and. mcomp == 33 .and. mphas == 11) then
        option%imode = RICHARDS_MODE
        option%mode = 'richards'
        option%nphase = 1; option%ndof = 2
        if (option%nspec > 1) then
          option%ndof = option%nspec +1
        endif
      endif
    endif
  endif
  if (option%imode == NULL_MODE) then
    call printErrMsg(option,"No mode specified")
  endif       

end subroutine setMode

! ************************************************************************** !
!
! assignMaterialPropToRegions: Assigns material properties to 
!                                    associated regions in the model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine assignMaterialPropToRegions(realization)

  use Realization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Fileio_module

  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: por0_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_pow_p(:)
  PetscReal, pointer :: tor_loc_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, material_id
  PetscInt :: istart, iend
  PetscInt :: fid = 86, status
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  
  type(material_type), pointer :: material
  type(region_type), pointer :: region
  
  option => realization%option
  grid => realization%grid
  field => realization%field

  ! loop over all strata to determine if any are inactive or
  ! have associated cell by cell material ids
  strata => realization%strata%first
  do
    if (.not.associated(strata)) exit
    if (.not.strata%active .or. .not.associated(strata%region)) then
      if (.not.associated(field%imat)) then
        allocate(field%imat(grid%ngmax))
        field%imat = -999
      endif
      exit
    endif
    strata => strata%next
  enddo
  
  ! Read in cell by cell material ids if they exist
  strata => realization%strata%first
  if (associated(strata)) then
    if (.not.associated(strata%region) .and. strata%active) then
      call readMaterialsFromFile(realization,strata%material_name)
    endif
  endif
    
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(field%porosity0,por0_p,ierr)
  call VecGetArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr)
  call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr)
  call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr)
  call VecGetArrayF90(field%perm_pow,perm_pow_p,ierr)

  strata => realization%strata%first
  do
    if (.not.associated(strata)) exit
   
    if (strata%active) then
      region => strata%region
      material => strata%material
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
        if (associated(field%imat)) then
          ! if field%imat is allocated and the id > 0, the material id 
          ! supercedes the material pointer for the strata
          material_id = field%imat(ghosted_id)
          if (material_id > 0 .and. &
              material_id <= size(realization%material_array)) then
            material => realization%material_array(material_id)%ptr
          endif
          ! otherwide set the imat value to the stratas material
          if (material_id < -998) & ! prevent overwrite of cell already set to inactive
            field%imat(ghosted_id) = material%id
        endif
        if (associated(material)) then
          icap_loc_p(ghosted_id) = material%icap
          ithrm_loc_p(ghosted_id) = material%ithrm
          por0_p(local_id) = material%porosity
          tor_loc_p(ghosted_id) = material%tortuosity
          perm_xx_p(local_id) = material%permeability(1,1)
          perm_yy_p(local_id) = material%permeability(2,2)
          perm_zz_p(local_id) = material%permeability(3,3)
          perm_pow_p(local_id) = material%permeability_pwr
        endif
      enddo
    endif
    strata => strata%next
  enddo

  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity0,por0_p,ierr)
  call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr)
  call VecRestoreArrayF90(field%perm_pow,perm_pow_p,ierr)
  call VecRestoreArrayF90(field%tor_loc,tor_loc_p,ierr)

  call GridGlobalToLocal(realization%grid,field%porosity0, &
                         field%porosity_loc,ONEDOF)
  call GridGlobalToLocal(realization%grid,field%perm0_xx, &
                         field%perm_xx_loc,ONEDOF)  
  call GridGlobalToLocal(realization%grid,field%perm0_yy, &
                         field%perm_yy_loc,ONEDOF)  
  call GridGlobalToLocal(realization%grid,field%perm0_zz, &
                         field%perm_zz_loc,ONEDOF)   

  call GridLocalToLocal(realization%grid,field%icap_loc, &
                        field%icap_loc,ONEDOF)   
  call GridLocalToLocal(realization%grid,field%ithrm_loc, &
                        field%ithrm_loc,ONEDOF)   
  call GridLocalToLocal(realization%grid,field%tor_loc, &
                        field%tor_loc,ONEDOF)   
  
end subroutine assignMaterialPropToRegions

! ************************************************************************** !
!
! assignInitialConditions: Assigns initial conditions to model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine assignInitialConditions(realization)

  use Realization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  
  use MPHASE_module, only : pflow_mphase_setupini

  implicit none
  
  type(realization_type) :: realization
  
  PetscReal, pointer :: pressure_p(:)
  PetscReal, pointer :: temp_p(:)
  PetscReal, pointer :: sat_p(:)
  PetscReal, pointer :: conc_p(:)
  PetscReal, pointer :: xmol_p(:)
  
  PetscInt :: icell, iconn, count, jn1, jn2, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: initial_condition
    
  option => realization%option
  field => realization%field
  grid => realization%grid
    
  select case(option%imode)
    case(RICHARDS_LITE_MODE)
    case(RICHARDS_MODE)
    case(MPH_MODE)
      call pflow_mphase_setupini(realization)
  end select 

  ! assign initial conditions values to domain
  call VecGetArrayF90(field%xx,xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  xx_p = -1.d20
  
  initial_condition => realization%initial_conditions%first
  do
  
    if (.not.associated(initial_condition)) exit

    if (.not.associated(initial_condition%connection)) then
      do icell=1,initial_condition%region%num_cells
        local_id = initial_condition%region%cell_ids(icell)
        ghosted_id = realization%grid%nL2G(local_id)
        iend = local_id*option%ndof
        ibegin = iend-option%ndof+1
        if (associated(field%imat)) then
          if (field%imat(ghosted_id) <= 0) then
            xx_p(ibegin:iend) = 0.d0
            iphase_loc_p(ghosted_id) = 0
            cycle
          endif
        endif
        do idof = 1, option%ndof
          xx_p(ibegin+idof) = &
            initial_condition%condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
        enddo
        iphase_loc_p(ghosted_id)=initial_condition%condition%iphase
      enddo
    else
      do iconn=1,initial_condition%connection%num_connections
        local_id = initial_condition%connection%id_dn(iconn)
        ghosted_id = realization%grid%nL2G(local_id)
        iend = local_id*option%ndof
        ibegin = iend-option%ndof+1
        if (associated(field%imat)) then
          if (field%imat(ghosted_id) <= 0) then
            xx_p(ibegin:iend) = 0.d0
            iphase_loc_p(ghosted_id) = 0
            cycle
          endif
        endif
        xx_p(ibegin:iend) = &
          initial_condition%aux_real_var(1:option%ndof,iconn)
        iphase_loc_p(ghosted_id)=initial_condition%aux_int_var(1,iconn)
      enddo
    endif
    initial_condition => initial_condition%next
  enddo
  
  call VecRestoreArrayF90(field%xx,xx_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  ! update dependent vectors
  call GridGlobalToLocal(grid,field%xx,field%xx_loc,NDOF)  
  call VecCopy(field%xx, field%yy, ierr)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)  
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_old_loc,ONEDOF)

end subroutine assignInitialConditions

! ************************************************************************** !
!
! verifyCouplers: Verifies the connectivity of a coupler
! author: Glenn Hammond
! date: 1/8/08
!
! ************************************************************************** !
subroutine verifyCouplers(realization,coupler_list)

  use Realization_module
  use Option_module 
  use Coupler_module
  use Grid_module
  use Output_module

  implicit none

  type(realization_type) :: realization
  type(coupler_list_type) :: coupler_list

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: coupler
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: iconn, local_id
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  grid => realization%grid

  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    call VecZeroEntries(global_vec,ierr)
    call VecGetArrayF90(global_vec,vec_ptr,ierr) 
    do iconn = 1, coupler%connection%num_connections
      local_id = coupler%connection%id_dn(iconn)
      vec_ptr(local_id) = coupler%id
    enddo
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr) 
    dataset_name = trim(coupler%condition%name) // '_' // &
                   trim(coupler%region%name)
    filename = trim(dataset_name) // '.tec'
    call OutputVectorTecplot(filename,dataset_name,realization,global_vec)

    coupler => coupler%next
  enddo

  call VecDestroy(global_vec,ierr)

end subroutine verifyCouplers

! ************************************************************************** !
!
! readRegionFiles: Reads in grid cell ids stored in files
! author: Glenn Hammond
! date: 1/03/08
!
! ************************************************************************** !
subroutine readRegionFiles(realization)

  use Realization_module
  use Region_module
  use HDF5_module

  implicit none

  type(realization_type) :: realization
  
  type(region_type), pointer :: region
  
  region => realization%regions%first
  do 
    if (.not.associated(region)) exit
    if (len_trim(region%filename) > 1) then
      if (index(region%filename,'.h5') > 0) then
        call HDF5ReadRegionFromFile(realization,region,region%filename)
      else
        call RegionReadFromFile(region,region%filename)
      endif
    endif
    region => region%next
  enddo

end subroutine readRegionFiles

! ************************************************************************** !
!
! readMaterialsFromFile: Reads in grid cell materials
! author: Glenn Hammond
! date: 1/03/08
!
! ************************************************************************** !
subroutine readMaterialsFromFile(realization,filename)

  use Realization_module
  use Field_module
  use Grid_module
  use Option_module
  use Fileio_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) realization
  character(len=MAXWORDLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr

  field => realization%field
  grid => realization%grid
  option => realization%option

  if (index(filename,'.h5') > 0) then
    call HDF5ReadMaterialsFromFile(realization,filename)
  else
    call GridCreateNaturalToGhostedHash(grid,option)
    status = 0
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      string = "File: " // filename // " not found."
      call printErrMsg(option,string)
    endif
    do
      call fiReadFlotranString(fid,string,ierr)
      if (ierr /= 0) exit
      call fiReadInt(string,natural_id,ierr)
      call fiErrorMsg(option%myrank,'natural id','STRATA', ierr)
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call fiReadInt(string,material_id,ierr)
        call fiErrorMsg(option%myrank,'material id','STRATA', ierr)
        field%imat(ghosted_id) = material_id
      endif
    enddo
    call GridDestroyHashTable(grid)
  endif
  
end subroutine readMaterialsFromFile

end module Init_module
