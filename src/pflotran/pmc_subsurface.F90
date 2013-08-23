module PMC_Subsurface_class

  use PMC_Base_class
  use Realization_class

  implicit none

#include "definitions.h"
  
  private

  type, public, extends(pmc_base_type) :: pmc_subsurface_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMCSubsurfaceInit
    procedure, public :: GetAuxData => PMCSubsurfaceGetAuxData
    procedure, public :: SetAuxData => PMCSubsurfaceSetAuxData
  end type pmc_subsurface_type
  
  public :: PMCSubsurfaceCreate
  
contains

! ************************************************************************** !
!
! PMCSubsurfaceCreate: Allocates and initializes a new process_model_coupler 
!                      object.
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function PMCSubsurfaceCreate()

  implicit none
  
  class(pmc_subsurface_type), pointer :: PMCSubsurfaceCreate
  
  class(pmc_subsurface_type), pointer :: pmc

  print *, 'PMCSubsurface%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceCreate => pmc  
  
end function PMCSubsurfaceCreate

! ************************************************************************** !
!
! PMCSubsurfaceInit: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PMCSubsurfaceInit(this)

  implicit none
  
  class(pmc_subsurface_type) :: this
  
  print *, 'PMCSubsurface%Init()'
  
  call PMCBaseInit(this)
  this%name = 'PMCSubsurface'
  nullify(this%realization)

end subroutine PMCSubsurfaceInit

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/22/13
! ************************************************************************** !
subroutine PMCSubsurfaceGetAuxData(this)

  use Connection_module
  use Coupler_module
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
!  use Realization_Base_class
  use Realization_class
  use String_module
  use Water_EOS_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_subsurface_type) :: this
  
  class(realization_type), pointer     :: realization
  type (patch_type),pointer            :: patch
  type (grid_type),pointer             :: grid
  type (coupler_list_type), pointer    :: coupler_list
  type (coupler_type), pointer         :: coupler
  type (option_type), pointer          :: option
  type (field_type),pointer            :: field
  type (connection_set_type), pointer  :: cur_connection_set
  PetscBool                            :: coupler_found
  PetscInt                             :: iconn
  PetscReal                            :: den
  PetscReal                            :: dt
  PetscReal, pointer                   :: mflux_p(:)
  PetscReal, pointer                   :: hflux_p(:)
  PetscErrorCode                       :: ierr

  print *, 'PMCSubsurfaceGetAuxData()'
  
  dt = this%option%surf_subsurf_coupling_flow_dt

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

      if (this%sim_aux%subsurf_mflux_exchange_with_surf /= 0) then
        ! PETSc Vector to store relevant mass-flux data between
        ! surface-subsurface model exists

        patch      => pmc%realization%patch
        grid       => pmc%realization%discretization%grid
        field      => pmc%realization%field
        option     => pmc%realization%option

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

            call density(option%reference_temperature, option%reference_pressure, &
                         den)

            coupler_list => patch%source_sinks
            coupler => coupler_list%first
            do
              if (.not.associated(coupler)) exit

              ! FLOW
              if (associated(coupler%flow_aux_real_var)) then

                ! Find the BC from the list of BCs
                if(StringCompare(coupler%name,'from_surface_ss')) then
                  coupler_found = PETSC_TRUE
                  
                  call VecGetArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                      mflux_p,ierr)
                  do iconn = 1,coupler%connection_set%num_connections
                    coupler%flow_aux_real_var(ONE_INTEGER,iconn) = -mflux_p(iconn)/dt*den
                  enddo
                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                          mflux_p,ierr)

                  call VecSet(pmc%sim_aux%surf_mflux_exchange_with_subsurf,0.d0,ierr)
                endif
              endif

              coupler => coupler%next
            enddo

          case (TH_MODE)
            this%option%io_buffer='Extend PMCSubsurfaceGetAuxData for TH'
            call printErrMsg(this%option)
        end select

      endif
    end select
  endif

end subroutine PMCSubsurfaceGetAuxData

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 08/21/13
! ************************************************************************** !
subroutine PMCSubsurfaceSetAuxData(this)

  use Grid_module
  use String_module
  use Realization_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Connection_module
  use Realization_Base_class

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_subsurface_type) :: this
  
  class(realization_type), pointer     :: realization
  type (patch_type),pointer            :: patch
  type (grid_type),pointer             :: grid
  type (coupler_list_type), pointer    :: coupler_list
  type (coupler_type), pointer         :: coupler
  type (option_type), pointer          :: option
  type (field_type),pointer            :: field
  type (connection_set_type), pointer  :: cur_connection_set
  PetscInt                             :: local_id
  PetscInt                             :: ghosted_id
  PetscInt                             :: iconn
  PetscInt                             :: istart
  PetscInt                             :: iend
  PetscReal, pointer                   :: xx_loc_p(:)
  PetscReal, pointer                   :: pres_top_bc_p(:)
  PetscReal, pointer                   :: temp_top_bc_p(:)
  PetscErrorCode                       :: ierr

  print *, 'PMCSubsurfaceSetAuxData()'

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

        if (this%sim_aux%subsurf_pres_top_bc/=0) then
          ! PETSc Vector to store relevant subsurface-flow data for
          ! surface-flow model exists

          patch      => pmc%realization%patch
          grid       => pmc%realization%discretization%grid
          field      => pmc%realization%field
          option     => pmc%realization%option

          coupler_list => pmc%realization%patch%source_sinks
          coupler => coupler_list%first
          do
            if (.not.associated(coupler)) exit

            ! FLOW
            if (associated(coupler%flow_aux_real_var)) then
              cur_connection_set => coupler%connection_set

              ! Find the BC from the list of BCs
              if (StringCompare(coupler%name,'from_surface_ss')) then

                select case(this%option%iflowmode)
                  case (RICHARDS_MODE)
                    call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc,pres_top_bc_p,ierr)
                    do iconn = 1,cur_connection_set%num_connections
                      local_id = cur_connection_set%id_dn(iconn)
                      ghosted_id = grid%nL2G(local_id)
                      pres_top_bc_p(iconn) = xx_loc_p(ghosted_id)
                    enddo
                    call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc,pres_top_bc_p,ierr)

                  case (TH_MODE)
                    call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc,pres_top_bc_p,ierr)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc,temp_top_bc_p,ierr)
                    do iconn = 1,cur_connection_set%num_connections
                      local_id = cur_connection_set%id_dn(iconn)
                      ghosted_id = grid%nL2G(local_id)
                      iend = ghosted_id*option%nflowdof
                      istart = iend-option%nflowdof+1

                      pres_top_bc_p(iconn) = xx_loc_p(istart)
                      temp_top_bc_p(iconn) = xx_loc_p(iend)
                    enddo
                    call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc,pres_top_bc_p,ierr)
                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc,temp_top_bc_p,ierr)
                  case default
                    call printErrMsg(this%option, &
                      'PMCSubsurfaceSetAuxData not supported for this MODE')
                end select
              endif

            endif

            coupler => coupler%next
          enddo

        endif
    end select

  endif

end subroutine PMCSubsurfaceSetAuxData


! ************************************************************************** !
!
! PMCSubsurfaceFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCSubsurfaceFinalizeRun(this)

  use Option_module
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
  call printMsg(this%option,'PMCSubsurface%FinalizeRun()')
  
  nullify(this%realization)
  
end subroutine PMCSubsurfaceFinalizeRun

! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
recursive subroutine Destroy(this)

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none
  
  class(pmc_subsurface_type) :: this
  
  call printMsg(this%option,'PMCSubsurface%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
end subroutine Destroy
  
end module PMC_Subsurface_class
