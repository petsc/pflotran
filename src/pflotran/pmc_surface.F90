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
  this%Synchronize1 => PMCSurfaceSynchronize1
  this%Synchronize2 => PMCSurfaceSynchronize2
  this%Synchronize3 => PMCSurfaceSynchronize3
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
  PetscErrorCode :: ierr
  
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call printVerboseMsg(this%option)
  
  if (associated(this%Synchronize1)) then
    !geh: PGI requires cast to base
    !call this%Synchronize1()
    pmc_base => this%CastToBase()
    !call pmc_base%Synchronize1()
  endif

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
    end select
    call this%timestepper%SetTargetTime(sync_time,this%option, &
                                        local_stop_flag,plot_flag, &
                                        transient_plot_flag)

    this%option%surf_flow_dt = this%timestepper%dt
    if (associated(this%Synchronize2)) then
      !geh: PGI requires cast to base
      pmc_base => this%CastToBase()
      !call pmc_base%Synchronize2()
    endif

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

  enddo
  
  this%option%surf_flow_time = this%timestepper%target_time

  if (associated(this%Synchronize3)) then
    !call this%Synchronize3()
  endif

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
!! date: 07/08/13
! ************************************************************************** !
subroutine PMCSurfaceSynchronize1(this)

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_base_type), pointer :: this
  PetscErrorCode :: ierr

  PetscReal :: dt

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED) then
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

  if (this%option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
    dt = this%option%surf_subsurf_coupling_flow_dt
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowUpdateSurfState(pmc%subsurf_realization, &
                                             pmc%surf_realization, dt)
          case (TH_MODE)
            call printErrMsg(this%option, &
            'Added code in PMCSurfaceSynchronize1 for SEQ_COUPLED_NEW in TH')
        end select
    end select
  endif

end subroutine PMCSurfaceSynchronize1

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/08/13
! ************************************************************************** !
subroutine PMCSurfaceSynchronize2(this)

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_base_type), pointer :: this
  PetscErrorCode :: ierr

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowSurf2SubsurfFlux(pmc%subsurf_realization, &
                                             pmc%surf_realization)
          case (TH_MODE)
            call SurfaceTHSurf2SubsurfFlux(pmc%subsurf_realization, &
                                           pmc%surf_realization)
        end select
    end select
  endif

end subroutine PMCSurfaceSynchronize2

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/08/13
! ************************************************************************** !
subroutine PMCSurfaceSynchronize3(this)

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_base_type), pointer :: this
  PetscReal :: dt
  PetscErrorCode :: ierr

  dt = this%option%surf_subsurf_coupling_flow_dt

  if(this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowUpdateSubsurfSS(pmc%subsurf_realization, &
                                            pmc%surf_realization, &
                                            dt)
          case (TH_MODE)
            call SurfaceTHUpdateSubsurfSS(pmc%subsurf_realization, &
                                            pmc%surf_realization,dt)
          end select
    end select
  endif
  
  if (this%option%subsurf_surf_coupling == SEQ_COUPLED_NEW) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call SurfaceFlowUpdateSubsurfBC(pmc%subsurf_realization, &
                                            pmc%surf_realization)
          case (TH_MODE)
            call printErrMsg(this%option, &
              'Added code in PMCSurfaceSynchronize3 for SEQ_COUPLED_NEW in TH mode')
          end select
    end select
  endif

end subroutine PMCSurfaceSynchronize3


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
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  PetscErrorCode :: ierr

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
            call SurfaceTHUpdateSurfBC(pmc%subsurf_realization, &
                                           pmc%surf_realization)
            this%option%io_buffer='Extend PMCSurfaceGetAuxData for TH'
            call printErrMsg(this%option)
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

  use Surface_Flow_module
  use Surface_TH_module
  use Option_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_surface_type) :: this

  PetscReal :: dt
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
            call SurfaceTHUpdateSubsurfSS(pmc%subsurf_realization, &
                                            pmc%surf_realization,dt)
            this%option%io_buffer='Extend PMCSurfaceGetAuxData for TH'
            call printErrMsg(this%option)
        end select
    end select
  endif

end subroutine PMCSurfaceSetAuxData

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
