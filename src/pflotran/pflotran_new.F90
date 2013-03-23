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
! (509) 375-3875
! glenn.hammond@pnnl.gov
! Richland, WA

!=======================================================================
program pflotran
  
  use Simulation_module
  use Option_module
  use Input_module
  use Logging_module
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscLogDouble :: timex_wall(4)

  PetscBool :: pflotranin_option_found
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  type(simulation_type), pointer :: simulation
  type(option_type), pointer :: option
  
!  nullify(stochastic)
  option => OptionCreate()
  option%fid_out = OUT_UNIT

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
  option%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename, &
                                 pflotranin_option_found,option)
  
  ! begin single simulation section

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    
  call LoggingCreate()

  simulation => SimulationCreate(option)

  call OptionCheckCommandLine(option)

  call PetscGetTime(timex_wall(1), ierr)
  option%start_time = timex_wall(1)

  call Simulation%Initialize()
  call Simulation%InitializeRun()
  call Simulation%ExecuteRun()
  call Simulation%FinalizeRun()

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
    
  ! end single simulation section
  
  call PetscOptionsSetValue('-options_left','no',ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue('-objects_left','yes',ierr)

  call OptionDestroy(option)
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)
  call exit(86)

end program pflotran
