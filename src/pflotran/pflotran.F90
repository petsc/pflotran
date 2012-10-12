!=======================================================================
! PFLOTRAN v2.0 LA-CC-09-047
!=======================================================================

!Copyright 2009. Los Alamos National Security, LLC. This material was produced under U.S. 
!Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated 
!by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has 
!rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS 
!NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE 
!USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software 
!should be clearly marked, so as not to confuse it with the version available from LANL.
!Additionally, this library is free software; you can redistribute it and/or modify it under the 
!terms of the GNU Lesser General Public License as published by the Free Software Foundation; 
!either version 2.1 of the License, or (at your option) any later version. Accordingly, this 
!library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
!the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser 
!General Public License for more details.

! Send all bug reports/questions/comments to:
!
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-16, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM

! or

! Glenn E. Hammond
! Pacific Northwest National Laboratory
! Energy and Environment Directorate
! MSIN K9-36
! (509) 395-3895
! glenn.hammond@pnl.gov
! Richland, WA

!=======================================================================
  program pflotran
  
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Input_module
  use Init_module
  use Logging_module
  use Stochastic_module
  use Stochastic_Aux_module
  use Regression_module
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscLogDouble :: timex_wall(4)

  PetscBool :: truth
  PetscBool :: option_found  
  PetscBool :: single_inputfile
  PetscInt :: i
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)
  type(stochastic_type), pointer :: stochastic
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  nullify(stochastic)
  option => OptionCreate()
  option%fid_out = OUT_UNIT
  single_inputfile = PETSC_TRUE

  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,option%global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,option%global_group,ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group
  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename,option_found,option)

  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)

  string = '-v'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) option%verbosity = 1

  string = '-multisimulation'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) then
    single_inputfile = PETSC_FALSE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) stochastic => StochasticCreate()

  call InitReadStochasticCardFromInput(stochastic,option)

  if (associated(stochastic)) then
    call StochasticInit(stochastic,option)
    call StochasticRun(stochastic,option)
  else

    if (single_inputfile) then
      PETSC_COMM_WORLD = MPI_COMM_WORLD
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    else
      call InitReadInputFilenames(option,filenames)
      temp_int = size(filenames) 
      call SimulationCreateProcessorGroups(option,temp_int)
      option%input_filename = filenames(option%mygroup_id)
      i = index(option%input_filename,'.',PETSC_TRUE)
      if (i > 1) then
        i = i-1
      else
        ! for some reason len_trim doesn't work on MS Visual Studio in 
        ! this location
        i = len(trim(option%input_filename)) 
      endif
      option%global_prefix = option%input_filename(1:i)
      write(string,*) option%mygroup_id
      option%group_prefix = 'G' // trim(adjustl(string))
    endif
 
    if (option%verbosity > 0) then 
      call PetscLogBegin(ierr)
      string = '-log_summary'
      call PetscOptionsInsertString(string, ierr)
    endif
    call LoggingCreate()

    simulation => SimulationCreate(option)
    realization => simulation%realization

    call OptionCheckCommandLine(option)

    call PetscGetTime(timex_wall(1), ierr)
    option%start_time = timex_wall(1)

    call Init(simulation)

#ifdef SURFACE_FLOW
    call StepperRun(simulation%realization,simulation%surf_realization, &
                    simulation%flow_stepper, &
                    simulation%tran_stepper, &
                    simulation%surf_flow_stepper)
#else
    call StepperRun(simulation%realization,simulation%flow_stepper, &
                    simulation%tran_stepper)
#endif

    call RegressionOutput(simulation%regression,simulation%realization, &
                          simulation%flow_stepper,simulation%tran_stepper)

  ! Clean things up.
    call SimulationDestroy(simulation)

  ! Final Time
    call PetscGetTime(timex_wall(2), ierr)
    
    if (option%myrank == option%io_rank) then

      if (option%print_to_screen) then
        write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
          (timex_wall(2)-timex_wall(1))/3600.d0
      endif
      if (option%print_to_file) then
        write(option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
          (timex_wall(2)-timex_wall(1))/3600.d0
      endif
    endif

    if (option%myrank == option%io_rank .and. option%print_to_file) then
      close(option%fid_out)
    endif

    call LoggingDestroy()
    
  endif
  
  call PetscOptionsSetValue('-options_left','no',ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue('-objects_left','yes',ierr)

  call OptionDestroy(option)
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)

  end program pflotran
