#ifdef SURFACE_FLOW
module Surf_Subsurf_Factory_module

  use Surf_Subsurf_Simulation_class

  implicit none

  private

#include "definitions.h"

  public :: SurfSubsurfaceInitialize

contains
! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceInitialize(simulation_base,option)

  use Option_module
  use Input_module
  use Timestepper_Base_class
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(surfsubsurface_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SurfSubsurfaceSimulationCreate(option)
  call SurfSubsurfaceInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine SurfSubsurfaceInitialize

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceInitializePostPETSc(simulation, option)

  use Simulation_module
  use Surface_Simulation_class
  use Subsurface_Simulation_class
  use Surface_Factory_module
  use Subsurface_Factory_module
  use Option_module
  use Init_module
  use Surface_Flow_module
  use Surface_TH_module
  use Simulation_Aux_module
  use PMC_Base_class
  
  implicit none
  
  class(surfsubsurface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  type(surface_simulation_type) :: surf_simulation
  type(subsurface_simulation_type) :: subsurf_simulation
  type(simulation_type), pointer :: simulation_old
  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscInt :: init_status
  
  ! process command line arguments specific to subsurface
  !call SurfSubsurfInitCommandLineSettings(option)
  
  allocate(simulation_old)
  simulation_old => SimulationCreate(option)
  call Init(simulation_old)

  call HijackSimulation(simulation_old,subsurf_simulation)
  call SubsurfaceJumpStart(subsurf_simulation)

   simulation%realization => simulation_old%realization
   simulation%flow_process_model_coupler => &
        subsurf_simulation%flow_process_model_coupler
   simulation%rt_process_model_coupler => &
        subsurf_simulation%rt_process_model_coupler
   simulation%regression => simulation_old%regression

  if(option%nsurfflowdof>0) then
    ! Both, Surface-Subsurface flow active
    call HijackSurfaceSimulation(simulation_old,surf_simulation)
    call SurfaceJumpStart(surf_simulation)

    simulation%process_model_coupler_list => &
      surf_simulation%process_model_coupler_list
    surf_simulation%process_model_coupler_list%next => &
      subsurf_simulation%process_model_coupler_list
    surf_simulation%surf_flow_process_model_coupler%subsurf_realization => &
      simulation_old%realization
    simulation%flow_process_model_coupler%realization => &
      simulation_old%realization

    simulation%surf_realization => simulation_old%surf_realization
    simulation%surf_flow_process_model_coupler => &
         surf_simulation%surf_flow_process_model_coupler

    nullify(surf_simulation%process_model_coupler_list)
  
    if (option%subsurf_surf_coupling == SEQ_COUPLED .or. &
        option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
       select case(option%iflowmode)
         case (RICHARDS_MODE)
            call SurfaceFlowGetSubsurfProp(simulation%realization, &
                 simulation%surf_realization)
         case (TH_MODE)
            call SurfaceTHGetSubsurfProp(simulation%realization, &
                 simulation%surf_realization)
         end select
      end if
   else
      ! Only subsurface flow active
      simulation%process_model_coupler_list => &
           subsurf_simulation%process_model_coupler_list
      ! call printErrMsg(option,'Only subsurface-flow is active. ' // &
      !        'Check inputfile or switch -simulation_mode subsurface')
      nullify(simulation%surf_realization)
      nullify(simulation%surf_flow_process_model_coupler)
   endif

   nullify(subsurf_simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors
  !call SimAuxCreateSubSurfVecs(simulation%sim_aux,simulation%realization, &
  !                                 option)
  !call SimAuxCreateSurfVecs(simulation%sim_aux,simulation%surf_realization, &
  !                              option)
  call SimAuxCreateSubSurfVecs(simulation%sim_aux,simulation_old%realization, &
                                   option)
  call SimAuxCreateSurfVecs(simulation%sim_aux,simulation_old%surf_realization, &
                                option)

  ! sim_aux: Create PETSc VectorScatters
  !call SimAuxCreateVecScatters(sim%sim_aux, &
  !                                 simulation%realization, &
  !                                 simulation%surf_realization)
  call SimAuxCreateVecScatters(simulation%sim_aux, &
                                   simulation_old%realization, &
                                   simulation_old%surf_realization)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if(associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if(associated(surf_simulation%surf_flow_process_model_coupler)) &
    surf_simulation%surf_flow_process_model_coupler%sim_aux => simulation%sim_aux

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%next)) then
    cur_process_model_coupler => cur_process_model_coupler%next
    call cur_process_model_coupler%SetAuxData()
  endif

  deallocate(simulation_old)

end subroutine SurfSubsurfaceInitializePostPETSc


! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/01/13
! ************************************************************************** !
subroutine SurfSubsurfInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
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
  
end subroutine SurfSubsurfInitCommandLineSettings

end module Surf_Subsurf_Factory_module

#endif
