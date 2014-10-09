module Geomechanics_Factory_module

  use Geomechanics_Simulation_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: GeomechanicsInitialize

contains

! ************************************************************************** !

subroutine GeomechanicsInitialize(simulation_base,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module
  use Input_Aux_module
  use Timestepper_Base_class
  use Simulation_Base_class

  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(geomechanics_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => GeomechanicsSimulationCreate(option)
  call GeomechanicsInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine GeomechanicsInitialize

! ************************************************************************** !

subroutine GeomechanicsInitializePostPETSc(simulation, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Init_module
  use Option_module
  use PMC_Base_class
  use PMC_Geomechanics_class
  use PFLOTRAN_Constants_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Realization_class
  use Geomechanics_Simulation_class
  use Simulation_module
  use Simulation_Aux_module
  use Subsurface_Simulation_class
  use Subsurface_Factory_module

  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(geomechanics_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  type(geomechanics_simulation_type) :: geomech_simulation
  type(subsurface_simulation_type) :: subsurf_simulation
  type(simulation_type), pointer :: simulation_old
  class(pmc_base_type), pointer :: cur_process_model_coupler
  type(gmdm_ptr_type), pointer                 :: dm_ptr

  allocate(simulation_old)
  simulation_old => SimulationCreate(option)
  call Init(simulation_old)

  call HijackSimulation(simulation_old, subsurf_simulation)
  call SubsurfaceJumpStart(subsurf_simulation)

  simulation%realization => simulation_old%realization
  simulation%flow_process_model_coupler => &
      subsurf_simulation%flow_process_model_coupler
  simulation%rt_process_model_coupler => &
      subsurf_simulation%rt_process_model_coupler
  simulation%regression => simulation_old%regression

  if (option%ngeomechdof > 0) then
    ! Both, Subsurface-flow and Geomechanics are active
    call HijackGeomechanicsSimulation(simulation_old, geomech_simulation)
    call GeomechanicsJumpStart(geomech_simulation)

    ! 1st PMC is subsurface-flow
    simulation%process_model_coupler_list => &
      subsurf_simulation%process_model_coupler_list

    ! 2st PMC is geomechanics
    subsurf_simulation%process_model_coupler_list%below => &
      geomech_simulation%process_model_coupler_list

    geomech_simulation%geomech_process_model_coupler%subsurf_realization => &
      simulation_old%realization
    simulation%flow_process_model_coupler%realization => &
      simulation_old%realization
    simulation%process_model_coupler_list%is_master = PETSC_TRUE
    geomech_simulation%process_model_coupler_list%is_master = PETSC_FALSE

    simulation%geomech_realization => simulation_old%geomech_realization
    simulation%geomech_process_model_coupler => &
         geomech_simulation%geomech_process_model_coupler

    nullify(geomech_simulation%process_model_coupler_list)

  else

    ! Only subsurface flow active
    simulation%process_model_coupler_list => &
         subsurf_simulation%process_model_coupler_list
    simulation%process_model_coupler_list%is_master = PETSC_TRUE
    ! call printErrMsg(option,'Only subsurface-flow is active. ' // &
    !        'Check inputfile or switch -simulation_mode subsurface')
    nullify(simulation%geomech_realization)
    nullify(simulation%geomech_process_model_coupler)

  endif

   nullify(subsurf_simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors and VectorScatters
  if(option%ngeomechdof > 0) then

    call GeomechCreateGeomechSubsurfVec(simulation_old%realization, &
                                        simulation_old%geomech_realization)
    call SimAuxCopySubsurfVec(simulation%sim_aux, simulation_old%realization%field%work)

    call GeomechCreateSubsurfStressStrainVec(simulation_old%realization, &
                                             simulation_old%geomech_realization)
    call SimAuxCopySubsurfGeomechVec(simulation%sim_aux, &
          simulation_old%geomech_realization%geomech_field%strain_subsurf)

    call GeomechRealizMapSubsurfGeomechGrid(simulation_old%realization, &
                                            simulation_old%geomech_realization, &
                                            option)

    dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
                simulation%geomech_realization%geomech_discretization, ONEDOF)

    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                              SUBSURF_TO_GEOMECHANICS)
    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                              GEOMECHANICS_TO_SUBSURF)
  endif

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if(associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if(option%ngeomechdof>0 .and. &
     associated(geomech_simulation%geomech_process_model_coupler)) &
    geomech_simulation%geomech_process_model_coupler%sim_aux => simulation%sim_aux

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%below)) then
    cur_process_model_coupler => cur_process_model_coupler%below
    call cur_process_model_coupler%GetAuxData()
    call cur_process_model_coupler%SetAuxData()
    select type(pmc => cur_process_model_coupler)
      class is(pmc_geomechanics_type)
        call GeomechStoreInitialPressTemp(pmc%geomech_realization)
    end select
  endif

  deallocate(simulation_old)

end subroutine GeomechanicsInitializePostPETSc

! ************************************************************************** !

subroutine HijackGeomechanicsSimulation(simulation_old,simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Simulation_module
  use Geomechanics_Realization_class
  use Option_module
  
  use PMC_Base_class
  use PMC_Geomechanics_class
  use Simulation_Base_class
  use PM_Geomechanics_Force_class
  use PM_Base_class
  use PM_Base_Pointer_module
  use Timestepper_Geomechanics_class

  implicit none
  
  type(simulation_type) :: simulation_old
  class(geomechanics_simulation_type) :: simulation
  
  class(pmc_geomechanics_type), pointer :: geomech_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(geomech_realization_type), pointer :: geomech_realization
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  geomech_realization => simulation_old%geomech_realization
  option => geomech_realization%option

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => geomech_realization%output_option
  simulation%option => geomech_realization%option
  nullify(cur_process_model)

  nullify(geomech_process_model_coupler)

  ! Create Surface-flow ProcessModel & ProcessModelCoupler
  if (option%ngeomechdof > 0) then
    cur_process_model => PMGeomechForceCreate()
    cur_process_model%option => geomech_realization%option
    cur_process_model%output_option => geomech_realization%output_option

    geomech_process_model_coupler => PMCGeomechanicsCreate()
    geomech_process_model_coupler%option => option
    geomech_process_model_coupler%pms => cur_process_model
    geomech_process_model_coupler%pm_ptr%ptr => cur_process_model
    call HijackTimestepper(simulation_old%geomech_timestepper, &
                           geomech_process_model_coupler%timestepper)
    nullify(cur_process_model)
  endif

  ! Add the ProcessModelCouplers in a list
  if (associated(geomech_process_model_coupler)) then
    simulation%process_model_coupler_list => geomech_process_model_coupler%CastToBase()
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoint_list => geomech_realization%waypoint_list
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pms
      do
        if (.not.associated(cur_process_model)) exit
        select type(cur_process_model)
          class is (pm_geomech_force_type)
            call cur_process_model%PMGeomechForceSetRealization(geomech_realization)
            call cur_process_model_coupler%SetTimestepper( &
                    geomech_process_model_coupler%timestepper)
            !geomech_process_model_coupler%timestepper%dt = option%flow_dt
        end select

        call cur_process_model%Init()
        select type(cur_process_model)
          class is (pm_geomech_force_type)
            select type(ts => cur_process_model_coupler%timestepper)
              class is (timestepper_geomechanics_type)
                call SNESSetFunction( &
                               ts%solver%snes, &
                               cur_process_model%residual_vec, &
                               PMResidual, &
                               cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
                call SNESSetJacobian( &
                               ts%solver%snes, &
                               ts%solver%J, &
                               ts%solver%Jpre, &
                               PMJacobian, &
                               cur_process_model_coupler%pm_ptr, &
                                     ierr);CHKERRQ(ierr)
            end select
        end select
        cur_process_model => cur_process_model%next
      enddo
      cur_process_model_coupler => cur_process_model_coupler%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo

  simulation%geomech_realization => geomech_realization
  simulation%geomech_process_model_coupler => geomech_process_model_coupler
  simulation%regression => simulation_old%regression
  geomech_process_model_coupler%geomech_realization => geomech_realization

end subroutine HijackGeomechanicsSimulation

! ************************************************************************** !

subroutine GeomechanicsJumpStart(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Geomechanics_Realization_class
  use Option_module
  use Timestepper_Geomechanics_class
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module
  use Condition_Control_module

  implicit none

  type(geomechanics_simulation_type) :: simulation

  class(geomech_realization_type), pointer :: geomch_realization
  class(timestepper_geomechanics_type), pointer :: master_timestepper
  class(timestepper_geomechanics_type), pointer :: geomech_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: geomech_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  geomch_realization => simulation%geomech_realization

  select type(ts => simulation%geomech_process_model_coupler%timestepper)
    class is(timestepper_geomechanics_type)
      geomech_timestepper => ts
  end select
  nullify(master_timestepper)

  option => geomch_realization%option
  output_option => geomch_realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", &
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  master_timestepper => geomech_timestepper

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  geomech_read = PETSC_FALSE
  failure = PETSC_FALSE
  

end subroutine GeomechanicsJumpStart

! ************************************************************************** !

subroutine HijackTimestepper(timestepper_old,timestepper_base)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Timestepper_Geomechanics_class
  use Timestepper_Base_class
  use Timestepper_module

  implicit none
  
  type(timestepper_type), pointer :: timestepper_old
  class(timestepper_base_type), pointer :: timestepper_base
  
  class(timestepper_geomechanics_type), pointer :: timestepper
  
  timestepper => TimestepperGeomechanicsCreate()
  
  timestepper%steps = timestepper_old%steps
  timestepper%num_constant_time_steps = timestepper_old%num_constant_time_steps

  timestepper%max_time_step = timestepper_old%max_time_step
  timestepper%max_time_step_cuts = timestepper_old%max_time_step_cuts
  timestepper%constant_time_step_threshold = timestepper_old%constant_time_step_threshold
  timestepper%cumulative_time_step_cuts = timestepper_old%cumulative_time_step_cuts
  timestepper%cumulative_solver_time = timestepper_old%cumulative_solver_time

  timestepper%start_time = timestepper_old%start_time
  timestepper%start_time_step = timestepper_old%start_time_step
  timestepper%time_step_tolerance = timestepper_old%time_step_tolerance
  timestepper%target_time = timestepper_old%target_time
  
  timestepper%prev_dt = timestepper_old%prev_dt
!  timestepper%dt = timestepper_old%dt
  timestepper%dt_min = timestepper_old%dt_min
  timestepper%dt_max = timestepper_old%dt_max
  timestepper%cfl_limiter = timestepper_old%cfl_limiter
  timestepper%cfl_limiter_ts = timestepper_old%cfl_limiter_ts
  
  timestepper%time_step_cut_flag = timestepper_old%time_step_cut_flag

  timestepper%init_to_steady_state = timestepper_old%init_to_steady_state
  timestepper%steady_state_rel_tol = timestepper_old%steady_state_rel_tol
  timestepper%run_as_steady_state = timestepper_old%run_as_steady_state

  timestepper%solver => timestepper_old%solver
  nullify(timestepper_old%solver)
  timestepper%cur_waypoint => timestepper_old%cur_waypoint
  nullify(timestepper_old%cur_waypoint)
  
  
!  stepper%prev_waypoint => timestepper_old%prev_waypoint
!  nullify(timestepper_old%prev_waypoint)
  
!  stepper%revert_dt = timestepper_old%revert_dt
!  stepper%num_contig_revert_due_to_sync = &
!  timestepper_old%num_contig_revert_due_to_sync

  timestepper_base => timestepper

end subroutine HijackTimestepper

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !

end module Geomechanics_Factory_module
