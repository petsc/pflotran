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
  use Realization_class
  use Timestepper_module
  use Option_module
  use Init_module
  use Stochastic_module
  use Stochastic_Aux_module
  use PFLOTRAN_Factory_module
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscInt :: init_status
  PetscErrorCode :: ierr
  type(stochastic_type), pointer :: stochastic
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(stepper_type), pointer :: master_stepper
  type(option_type), pointer :: option
  
  call InitializePFLOTRAN(option)
  
  select case(option%simulation_type)
    case(STOCHASTIC_SIM_TYPE)
      stochastic => StochasticCreate()
      call StochasticInit(stochastic,option)
      call StochasticRun(stochastic,option)
    case(SUBSURFACE_SIM_TYPE,MULTISIMULATION_SIM_TYPE)
      if (option%simulation_type == SUBSURFACE_SIM_TYPE) then
        PETSC_COMM_WORLD = MPI_COMM_WORLD
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      else
        call InitDivvyUpSimulations(option)
      endif

      simulation => SimulationCreate(option)
      realization => simulation%realization

      call Init(simulation)

#ifdef SURFACE_FLOW
      call TimestepperInitializeRun(simulation%realization, &
                                    simulation%surf_realization, &
                                    master_stepper, &
                                    simulation%flow_stepper, &
                                    simulation%tran_stepper, &
                                    simulation%surf_flow_stepper, &
                                    init_status)
      select case(init_status)
        case(TIMESTEPPER_INIT_PROCEED)
          call  TimestepperExecuteRun(simulation%realization, &
                                      simulation%surf_realization, &
                                      master_stepper, &
                                      simulation%flow_stepper, &
                                      simulation%tran_stepper, &
                                      simulation%surf_flow_stepper)
          call  TimestepperFinalizeRun(simulation%realization, &
                                       simulation%surf_realization, &
                                       master_stepper, &
                                       simulation%flow_stepper, &
                                       simulation%tran_stepper, &
                                       simulation%surf_flow_stepper)
        case(TIMESTEPPER_INIT_FAIL)
        case(TIMESTEPPER_INIT_DONE)
      end select
#else
      call TimestepperInitializeRun(simulation%realization, &
                                    master_stepper, &
                                    simulation%flow_stepper, &
                                    simulation%tran_stepper, &
                                    init_status)
      select case(init_status)
        case(TIMESTEPPER_INIT_PROCEED)
          call  TimestepperExecuteRun(simulation%realization, &
                                      master_stepper, &
                                      simulation%flow_stepper, &
                                      simulation%tran_stepper)
          call  TimestepperFinalizeRun(simulation%realization, &
                                       master_stepper, &
                                       simulation%flow_stepper, &
                                       simulation%tran_stepper)
        case(TIMESTEPPER_INIT_FAIL)
        case(TIMESTEPPER_INIT_DONE)
      end select
#endif
  end select

  call FinalizePFLOTRAN(simulation)
  call exit(86)  

end program pflotran
