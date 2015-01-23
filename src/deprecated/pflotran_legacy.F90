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
! Sandia National Laboratories
! Applied Systems Analysis & Research
! 413 Cherry Blossom Lp
! (509) 392-1715
! gehammo@sandia.gov
! Richland, WA 99352

!=======================================================================
program pflotran
  
  use Option_module
  use Input_Aux_module
  use Stochastic_module
  use Stochastic_Aux_module
  use Simulation_module
  use Timestepper_module  
  use PFLOTRAN_Factory_module
  use Logging_module
  use Geomechanics_Logging_module
  
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)
  type(simulation_type), pointer :: simulation
  type(timestepper_type), pointer :: master_timestepper
  type(stochastic_type), pointer :: stochastic
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  PetscInt :: init_status
  
  option => OptionCreate()
  call OptionInitMPI(option)
  call PFLOTRANInitializePrePETSc(option)

  !geh: in the refactored code, this portion has been relocated to
  !     subsurface_factory.F90.
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif
  
  select case(option%subsurface_simulation_type)
    case(STOCHASTIC_SIM_TYPE)
      stochastic => StochasticCreate()
      ! PETSc initialized within StochasticInit()
      call StochasticInit(stochastic,option)
      call StochasticRun(stochastic,option)
    case(SUBSURFACE_SIM_TYPE,MULTISIMULATION_SIM_TYPE)
      if (option%subsurface_simulation_type == MULTISIMULATION_SIM_TYPE) then
        call InputReadFilenames(option,filenames)
        call OptionDivvyUpSimulations(option,filenames)
        deallocate(filenames)
        nullify(filenames)
      endif
      call OptionInitPetsc(option)
      call LoggingCreate()
      call GeomechLoggingCreate()
      call PFLOTRANInitializePostPETSc(simulation,master_timestepper,option, &
                                       init_status)
      call PFLOTRANRun(simulation,master_timestepper,init_status)
      call PFLOTRANFinalize(simulation,option)
      call LoggingDestroy()
      call GeomechLoggingDestroy()
  end select

  call OptionFinalize(option)

end program pflotran
