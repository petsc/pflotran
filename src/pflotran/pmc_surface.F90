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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
function PMCSurfaceCreate()

  implicit none
  
  class(pmc_surface_type), pointer :: PMCSurfaceCreate
  
  class(pmc_surface_type), pointer :: pmc

  print *, 'PMCSurfaceCreate%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSurfaceCreate => pmc  
  
end function PMCSurfaceCreate

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
subroutine PMCSurfaceInit(this)

  implicit none
  
  class(pmc_surface_type) :: this
  
  print *, 'PMCSurfaceInit%Init()'
  
  call PMCBaseInit(this)
  nullify(this%surf_realization)
  nullify(this%subsurf_realization)
!  nullify(this%surf_timestepper)

end subroutine PMCSurfaceInit

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
recursive subroutine PMCSurfaceRunToTime(this,sync_time,stop_flag)

  use Timestepper_Base_class
  use Output_module, only : Output
  use Realization_class, only : realization_type
  use Process_Model_Base_class
  use Process_Model_Surface_Flow_class
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
                                        transient_plot_flag)

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
    !if (associated(this%Output)) then
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
    !endif

    if (this%is_master .and. &
        this%option%checkpoint_flag .and. &
        mod(this%timestepper%steps, &
        this%option%checkpoint_frequency) == 0) then
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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/21/13
! ************************************************************************** !
subroutine PMCSurfaceAccumulateAuxData(this)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/21/13
! ************************************************************************** !
subroutine PMCSurfaceGetAuxData(this)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/21/13
! ************************************************************************** !
subroutine PMCSurfaceSetAuxData(this)

  use Grid_module
  use Option_module
  use Patch_module
  use Surface_Global_Aux_module
  use Surface_Flow_module
  use Surface_TH_module
  use Surface_TH_Aux_module
  use Surface_Realization_class

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  type(grid_type), pointer :: surf_grid
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars(:)
  type(patch_type), pointer :: surf_patch
  class(surface_realization_type), pointer :: surf_realization

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iend
  PetscInt :: istart

  PetscReal :: dt
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: surf_head_p(:)
  PetscReal, pointer :: surf_temp_p(:)
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
            surf_global_aux_vars => surf_patch%surf_aux%SurfaceGlobal%aux_vars
            surf_aux_vars => surf_patch%surf_aux%SurfaceTH%aux_vars

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                xx_loc_p,ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_head, surf_head_p, ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_temp, surf_temp_p, ierr)

            do ghosted_id = 1, surf_grid%ngmax
              local_id = surf_grid%nG2L(ghosted_id)
              if (local_id < 1) cycle
              iend = local_id*this%option%nflowdof
              istart = iend - this%option%nflowdof+1
              if (xx_loc_p(istart) < 1.d-15) then
                surf_head_p(local_id) = 0.d0
                surf_temp_p(local_id) = 0.d0
              else
                surf_head_p(local_id) = xx_loc_p(istart)
                surf_temp_p(local_id) = surf_global_aux_vars(ghosted_id)%temp(1)
              endif
            enddo

            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                    xx_loc_p,ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_head, surf_head_p,ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_temp, surf_temp_p,ierr)

        end select
    end select
  endif

end subroutine PMCSurfaceSetAuxData

! ************************************************************************** !
!> This routine is called to set values in sim_aux PETSc vectors after restart
!! checkpoint files is read.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 09/23/13
! ************************************************************************** !
subroutine PMCSurfaceGetAuxDataAfterRestart(this)

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module
  use Water_EOS_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: count
  PetscReal, pointer      :: hw_p(:)           ! head [m]
  PetscReal, pointer      :: surfpress_p(:)
  PetscReal :: den
  PetscErrorCode :: ierr

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

            call density(this%option%reference_temperature,this%option%reference_pressure,den)

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx, hw_p, ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p, ierr)
            count = 0
            do ghosted_id = 1, pmc%surf_realization%discretization%grid%ngmax

              local_id = pmc%surf_realization%discretization%grid%nG2L(ghosted_id)
              if(local_id <= 0) cycle

              count = count + 1
              surfpress_p(count) = hw_p(ghosted_id)*den*abs(this%option%gravity(3)) + &
                                   this%option%reference_pressure
            enddo
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx, hw_p, ierr)
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
            call SurfaceTHUpdateSurfBC(pmc%subsurf_realization, &
                                           pmc%surf_realization)
            this%option%io_buffer='Extend PMCSurfaceGetAuxData for TH'
            call printErrMsg(this%option)
        end select
    end select
  endif

end subroutine PMCSurfaceGetAuxDataAfterRestart

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
recursive subroutine PMCSurfaceFinalizeRun(this)

  use Option_module
  
  implicit none
  
  class(pmc_surface_type), target :: this
  
  call printMsg(this%option,'PMCSurface%FinalizeRun()')
  
  nullify(this%surf_realization)
!  nullify(this%surf_timestepper)
  
end subroutine PMCSurfaceFinalizeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
recursive subroutine Destroy(this)

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
