#ifdef SURFACE_FLOW
module PMC_Surface_class

  use PMC_Base_class
  use Realization_class
  use Surface_Realization_class
  use Timestepper_Surface_class

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  private

  type, public, extends(pmc_base_type) :: pmc_surface_type
    class(realization_type), pointer :: subsurf_realization
    class(surface_realization_type), pointer :: surf_realization
  contains
    procedure, public :: Init => PMCSurfaceInit
    procedure, public :: RunToTime => PMCSurfaceRunToTime
    procedure, public :: AccumulateAuxData => PMCSurfaceAccumulateAuxData
    procedure, public :: GetAuxData => PMCSurfaceGetAuxData
    procedure, public :: SetAuxData => PMCSurfaceSetAuxData
    procedure, public :: PMCSurfaceGetAuxDataAfterRestart
  end type pmc_surface_type

  public :: PMCSurfaceCreate

contains

! ************************************************************************** !

function PMCSurfaceCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(pmc_surface_type), pointer :: PMCSurfaceCreate
  
  class(pmc_surface_type), pointer :: pmc

  print *, 'PMCSurfaceCreate%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSurfaceCreate => pmc  
  
end function PMCSurfaceCreate

! ************************************************************************** !

subroutine PMCSurfaceInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(pmc_surface_type) :: this
  
  print *, 'PMCSurfaceInit%Init()'
  
  call PMCBaseInit(this)
  nullify(this%surf_realization)
  nullify(this%subsurf_realization)
!  nullify(this%surf_timestepper)

end subroutine PMCSurfaceInit

! ************************************************************************** !

recursive subroutine PMCSurfaceRunToTime(this,sync_time,stop_flag)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Timestepper_Base_class
  use Output_module, only : Output
  use Realization_class, only : realization_type
  use PM_Base_class
  use PM_Surface_Flow_class
  use Option_module
  use Surface_Flow_module
  use Surface_TH_module
  use Output_Surface_module
  
  implicit none
#include "finclude/petscviewer.h"
  
  class(pmc_surface_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  class(pmc_base_type), pointer :: pmc_base
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  PetscBool :: checkpoint_flag
  class(pm_base_type), pointer :: cur_pm
  PetscReal :: dt_max_loc
  PetscReal :: dt_max_glb
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = 0
  do
    if (local_stop_flag > 0) exit ! end simulation
    if (this%timestepper%target_time >= sync_time) exit
    
    call SetOutputFlags(this)
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE
    checkpoint_flag = PETSC_FALSE
    
    cur_pm => this%pm_list

    select case(this%option%iflowmode)
      case (RICHARDS_MODE)
        call SurfaceFlowComputeMaxDt(this%surf_realization,dt_max_loc)
      case (TH_MODE)
        call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
    end select

    ! Find mininum allowable timestep across all processors
    call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION,&
                     MPI_MIN,this%option%mycomm,ierr)

    select type(timestepper => this%timestepper)
      class is(timestepper_surface_type)
        timestepper%dt_max_allowable = dt_max_glb
        timestepper%surf_subsurf_coupling_flow_dt = &
          this%option%surf_subsurf_coupling_flow_dt
    end select
    call this%timestepper%SetTargetTime(sync_time,this%option, &
                                        local_stop_flag,plot_flag, &
                                        transient_plot_flag,checkpoint_flag)

    this%option%surf_flow_dt = this%timestepper%dt

    ! Accumulate data needed by process-model
    call this%AccumulateAuxData()

    call this%timestepper%StepDT(this%pm_list,local_stop_flag)

    if (local_stop_flag > 1) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_pm => this%pm_list
    do
      if (.not.associated(cur_pm)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_pm%UpdateSolution()
      !! TODO(gb)
      !!!call this%timestepper%UpdateDT(cur_pm)
      cur_pm => cur_pm%next
    enddo

#if 0
    ! Run underlying process model couplers
    if (associated(this%below)) then
      call this%below%RunToTime(this%timestepper%target_time,local_stop_flag)
    endif
#endif
    
    ! only print output for process models of depth 0
    ! TODO(GB): Modify OutputSurface()
    if (associated(this%Output)) then
      if (this%timestepper%time_step_cut_flag) then
        plot_flag = PETSC_FALSE
      endif
      ! however, if we are using the modulus of the output_option%imod, we may
      ! still print
      if (mod(this%timestepper%steps, &
              this%pm_list% &
                output_option%periodic_output_ts_imod) == 0) then
        plot_flag = PETSC_TRUE
      endif
      if (plot_flag .or. mod(this%timestepper%steps, &
                             this%pm_list%output_option% &
                               periodic_tr_output_ts_imod) == 0) then
        transient_plot_flag = PETSC_TRUE
      endif
      !call this%Output(this%pm_list%realization_base,plot_flag, &
      !                 transient_plot_flag)
      call OutputSurface(this%surf_realization, this%subsurf_realization, &
                         plot_flag, transient_plot_flag)
    endif

    if (this%is_master) then
      if (.not.checkpoint_flag) then
        if (this%option%checkpoint_flag .and. this%option%checkpoint_frequency > 0) then
          if (mod(this%timestepper%steps,this%option%checkpoint_frequency) == 0) then
           checkpoint_flag = PETSC_TRUE
          endif
        endif
       endif
    else
      checkpoint_flag = PETSC_FALSE
    endif

    if (checkpoint_flag) then
      ! if checkpointing, need to sync all other PMCs.  Those "below" are
      ! already in sync, but not those "next".
      ! Set data needed by process-model
      call this%SetAuxData()
      ! Run neighboring process model couplers
      if (associated(this%next)) then
        call this%next%RunToTime(this%timestepper%target_time,local_stop_flag)
      endif
      call this%GetAuxData()
      call this%Checkpoint(viewer,this%timestepper%steps)
    endif

  enddo
  
  this%option%surf_flow_time = this%timestepper%target_time

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%next)) then
    call this%next%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)
  
end subroutine PMCSurfaceRunToTime

! ************************************************************************** !

subroutine PMCSurfaceAccumulateAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module

  implicit none

  class(pmc_surface_type) :: this

  PetscErrorCode :: ierr

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowSurf2SubsurfFlux(pmc%subsurf_realization, &
                                             pmc%surf_realization)
            call VecCopy(pmc%surf_realization%surf_field%exchange_subsurf_2_surf, &
                         pmc%sim_aux%surf_mflux_exchange_with_subsurf,ierr)
            call VecSet(pmc%surf_realization%surf_field%exchange_subsurf_2_surf,0.d0,ierr)
          case (TH_MODE)
            call SurfaceTHSurf2SubsurfFlux(pmc%subsurf_realization, &
                                           pmc%surf_realization)
            this%option%io_buffer='Extend PMCSurfaceAccumulateAuxData for TH'
            call printErrMsg(this%option)
        end select
    end select
  endif

end subroutine PMCSurfaceAccumulateAuxData

! ************************************************************************** !

subroutine PMCSurfaceGetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

  use Surface_Flow_module
  use Surface_TH_module
  use Surface_TH_module
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  PetscErrorCode :: ierr

  print *, 'PMCSurfaceGetAuxData()'
  if (this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowUpdateSurfBC(pmc%subsurf_realization, &
                                             pmc%surf_realization)
          case (TH_MODE)
            call SurfaceTHUpdateSurfBC(pmc%subsurf_realization, &
                                           pmc%surf_realization)
        end select
    end select
  endif

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)
            call SurfaceFlowUpdateSurfStateNew(pmc%surf_realization)
          case (TH_MODE)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 pmc%surf_realization%surf_field%temp_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               pmc%surf_realization%surf_field%temp_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)
            call SurfaceTHUpdateSurfStateNew(pmc%surf_realization)
        end select
    end select
  endif

end subroutine PMCSurfaceGetAuxData

! ************************************************************************** !

subroutine PMCSurfaceSetAuxData(this)
  ! 
  ! This routine extracts data from surface flow model and stores it sim-aux,
  ! which will be required by the subsurface flow model.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

  use Connection_module
  use Coupler_module
  use Grid_module
  use Option_module
  use Patch_module
  use Surface_Global_Aux_module
  use Surface_Flow_module
  use Surface_TH_module
  use Surface_TH_Aux_module
  use Surface_Realization_class
  use String_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  type(grid_type), pointer :: surf_grid
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(patch_type), pointer :: surf_patch
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  class(surface_realization_type), pointer :: surf_realization

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iend
  PetscInt :: istart
  PetscInt :: iconn

  PetscReal :: dt
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: surf_head_p(:)
  PetscReal, pointer :: surf_temp_p(:)
  PetscReal, pointer :: surf_hflux_p(:)
  PetscBool :: found
  PetscReal :: esrc
  PetscErrorCode :: ierr

  dt = this%option%surf_subsurf_coupling_flow_dt

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_mflux_exchange_with_subsurf, &
                                 pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_mflux_exchange_with_subsurf, &
                               pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)

          case (TH_MODE)
            call SurfaceTHUpdateSubsurfSS(pmc%subsurf_realization, &
                                            pmc%surf_realization,dt)
            this%option%io_buffer='Extend PMCSurfaceGetAuxData for TH'
            call printErrMsg(this%option)
        end select
    end select
  endif


  if(this%option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
    select type(pmc => this)
      class is(pmc_surface_type)

        select case(this%option%iflowmode)

          case (RICHARDS_MODE)
            call VecCopy(pmc%surf_realization%surf_field%flow_xx, &
                         pmc%sim_aux%surf_head, ierr)
          case (TH_MODE)

            surf_realization => pmc%surf_realization
            surf_patch => surf_realization%patch
            surf_grid => surf_patch%grid
            surf_global_auxvars => surf_patch%surf_aux%SurfaceGlobal%auxvars
            surf_auxvars => surf_patch%surf_aux%SurfaceTH%auxvars

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                xx_loc_p,ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_head, surf_head_p, ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_temp, surf_temp_p, ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                surf_hflux_p, ierr)

            do ghosted_id = 1, surf_grid%ngmax
              local_id = surf_grid%nG2L(ghosted_id)
              if (local_id < 1) cycle
              iend = local_id*this%option%nflowdof
              istart = iend - this%option%nflowdof+1
              if (xx_loc_p(istart) < 1.d-15) then
                surf_head_p(local_id) = 0.d0
                surf_temp_p(local_id) = this%option%reference_temperature
              else
                surf_head_p(local_id) = xx_loc_p(istart)
                surf_temp_p(local_id) = surf_global_auxvars(ghosted_id)%temp(1)
              endif
            enddo

            found = PETSC_FALSE
            source_sink => surf_patch%source_sinks%first
            do
              if (.not.associated(source_sink)) exit

              if (associated(source_sink%flow_aux_real_var)) then
                cur_connection_set => source_sink%connection_set

                if (StringCompare(source_sink%name,'atm_energy_ss')) then

                  do iconn = 1, cur_connection_set%num_connections

                    local_id = cur_connection_set%id_dn(iconn)
                    select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
                      case (ENERGY_RATE_SS)
                        esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)
                      case (HET_ENERGY_RATE_SS)
                        esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
                      case default
                        this%option%io_buffer = 'atm_energy_ss does not have '// &
                          'a temperature condition that is either a ' // &
                          ' ENERGY_RATE_SS or HET_ENERGY_RATE_SS'
                    end select

                    ! Only when no standing water is present, the atmospheric
                    ! energy flux is applied directly on subsurface domain.
                    if (surf_head_p(local_id) == 0.d0) then
                      surf_hflux_p(local_id) = esrc
                    else
                      surf_hflux_p(local_id) = 0.d0
                    endif

                  enddo

                  found = PETSC_TRUE

                endif ! StringCompare()
              endif ! associate()

              source_sink => source_sink%next
            enddo

            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                    xx_loc_p,ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_head, surf_head_p,ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_temp, surf_temp_p,ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                surf_hflux_p, ierr)

            if (.not.(found)) then
              this%option%io_buffer = 'atm_energy_ss not found in surface-flow model'
              call printErrMsg(this%option)
            endif
        end select
    end select
  endif

end subroutine PMCSurfaceSetAuxData

! ************************************************************************** !

subroutine PMCSurfaceGetAuxDataAfterRestart(this)
  ! 
  ! This routine is called to set values in sim_aux PETSc vectors after restart
  ! checkpoint files is read.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/23/13
  ! 

  use Surface_Flow_module
  use Surface_TH_Aux_module
  use Surface_TH_module
  use Option_module
  use EOS_Water_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: count
  PetscReal, pointer      :: xx_p(:)
  PetscReal, pointer      :: surfpress_p(:)
  PetscReal, pointer      :: surftemp_p(:)
  PetscInt :: istart, iend
  PetscReal :: den
  PetscErrorCode :: ierr
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)

  print *, 'PMCSurfaceGetAuxDataAfterRestart()'
  if (this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
          case (TH_MODE)
        end select
    end select
  endif

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)

            call EOSWaterdensity(this%option%reference_temperature,this%option%reference_pressure,den)

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p, ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p, ierr)
            count = 0
            do ghosted_id = 1, pmc%surf_realization%discretization%grid%ngmax

              local_id = pmc%surf_realization%discretization%grid%nG2L(ghosted_id)
              if(local_id <= 0) cycle

              count = count + 1
              surfpress_p(count) = xx_p(ghosted_id)*den*abs(this%option%gravity(3)) + &
                                   this%option%reference_pressure
            enddo
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p, ierr)
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p, ierr)

            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 INSERT_VALUES,SCATTER_REVERSE,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               INSERT_VALUES,SCATTER_REVERSE,ierr)

          case (TH_MODE)

            ! NOTE(GB:) This is strictly not correct since density should be
            ! computed based on surface-water temperature (not on
            ! reference-temperature). Presently, SurfaceCheckpointProcessModel()
            ! does not output surface-water temperature for TH-Mode and the
            ! subroutine needs to be modified in future.
            call EOSWaterdensity(this%option%reference_temperature,this%option%reference_pressure,den)

            surf_auxvars => pmc%surf_realization%patch%surf_aux%SurfaceTH%auxvars

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p, ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p, ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%temp_subsurf, surftemp_p, ierr)

            count = 0
            do ghosted_id = 1, pmc%surf_realization%discretization%grid%ngmax

              local_id = pmc%surf_realization%discretization%grid%nG2L(ghosted_id)
              if(local_id <= 0) cycle

              count = count + 1
              iend = ghosted_id*this%option%nflowdof
              istart = iend - this%option%nflowdof+1
              surfpress_p(count) = xx_p(istart)*den*abs(this%option%gravity(3)) + &
                                   this%option%reference_pressure
              surftemp_p = xx_p(iend)/xx_p(istart)/den/ &
                      surf_auxvars(ghosted_id)%Cwi - 273.15d0
            enddo
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p, ierr)
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p, ierr)
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%temp_subsurf, surftemp_p, ierr)

            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 pmc%surf_realization%surf_field%temp_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               pmc%surf_realization%surf_field%temp_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD,ierr)
        end select
    end select
  endif

end subroutine PMCSurfaceGetAuxDataAfterRestart

! ************************************************************************** !

recursive subroutine PMCSurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_surface_type), target :: this
  
  call printMsg(this%option,'PMCSurface%FinalizeRun()')
  
  nullify(this%surf_realization)
!  nullify(this%surf_timestepper)
  
end subroutine PMCSurfaceFinalizeRun

! ************************************************************************** !

recursive subroutine Destroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none
  
  class(pmc_surface_type) :: this
  
  call printMsg(this%option,'PMCSurface%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
end subroutine Destroy

end module PMC_Surface_class
#endif
