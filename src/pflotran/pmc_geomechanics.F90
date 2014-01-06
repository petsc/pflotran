#ifdef GEOMECH

module PMC_Geomechanics_class

  use PMC_Base_class
  use Realization_class
#ifdef PROCESS_MODEL
  use Geomechanics_Realization_class
#else
  use Geomechanics_Realization_module
#endif
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  private

  type, public, extends(pmc_base_type) :: pmc_geomechanics_type
    class(realization_type), pointer :: subsurf_realization
    class(geomech_realization_type), pointer :: geomech_realization
  contains
    procedure, public :: Init => PMCGeomechanicsInit
    procedure, public :: RunToTime => PMCGeomechanicsRunToTime
    procedure, public :: GetAuxData => PMCGeomechanicsGetAuxData
    procedure, public :: SetAuxData => PMCGeomechanicsSetAuxData
  end type pmc_geomechanics_type

  public :: PMCGeomechanicsCreate

contains

! ************************************************************************** !
!> This routine allocates and initializes a new object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !
function PMCGeomechanicsCreate()

  implicit none
  
  class(pmc_geomechanics_type), pointer :: PMCGeomechanicsCreate
  
  class(pmc_geomechanics_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCGeomechanicsCreate%Create()'
#endif

  allocate(pmc)
  call pmc%Init()
  
  PMCGeomechanicsCreate => pmc  
  
end function PMCGeomechanicsCreate

! ************************************************************************** !
!> This routine initializes a new process model coupler object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !
subroutine PMCGeomechanicsInit(this)

  implicit none
  
  class(pmc_geomechanics_type) :: this
  
#ifdef DEBUG
  print *, 'PMCGeomechanics%Init()'
#endif

  call PMCBaseInit(this)
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsInit

! ************************************************************************** !
!> This routine runs the geomechanics simulation.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !
subroutine PMCGeomechanicsRunToTime(this,sync_time,stop_flag)

  use Timestepper_Base_class
  use Option_module
  use Process_Model_Base_class

  implicit none

  class(pmc_geomechanics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag

  PetscInt :: local_stop_flag
  class(pm_base_type), pointer :: cur_pm

  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = 0

  call this%timestepper%StepDT(this%pm_list,local_stop_flag)

  ! Have to loop over all process models coupled in this object and update
  ! the time step size.  Still need code to force all process models to
  ! use the same time step size if tightly or iteratively coupled.
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    ! Geomechanics PM does not have an associate time 
    !call this%timestepper%UpdateDT(cur_pm)
    cur_pm => cur_pm%next
  enddo

  ! Run underlying process model couplers
  if (associated(this%below)) then
    ! Set data needed by process-model
    call this%SetAuxData()
    call this%below%RunToTime(this%timestepper%target_time,local_stop_flag)
  endif

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%next)) then
    call this%next%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)

end subroutine PMCGeomechanicsRunToTime

! ************************************************************************** !
!> This routine updates data in simulation_aux that is required by other
!! process models.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !
subroutine PMCGeomechanicsSetAuxData(this)

  use Option_module
  use Grid_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(pmc_geomechanics_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(grid_type), pointer :: grid
  PetscInt :: local_id
  PetscScalar, pointer :: por0_p(:)
  PetscScalar, pointer :: por_p(:)
  PetscScalar, pointer :: strain_p(:)
  PetscErrorCode :: ierr
  PetscReal :: trace_epsilon
  PetscReal :: por_new

  ! If at initialization stage, do nothing
  if (this%timestepper%steps == 0) return

  select type(pmc => this)
    class is(pmc_geomechanics_type)
      if (this%option%geomech_subsurf_coupling == GEOMECH_TWO_WAY_COUPLED) then

        grid => pmc%subsurf_realization%patch%grid

        ! Save strain dataset in sim_aux%subsurf_strain
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%strain, &
                             pmc%sim_aux%subsurf_strain, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%strain, &
                           pmc%sim_aux%subsurf_strain, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
                             
        ! Save stress dataset in sim_aux%subsurf_stress
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%stress, &
                             pmc%sim_aux%subsurf_stress, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%stress, &
                           pmc%sim_aux%subsurf_stress, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! Update porosity dataset in sim_aux%subsurf_por
        call VecGetArrayF90(pmc%sim_aux%subsurf_por0, por0_p, ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_por, por_p, ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_strain, strain_p, ierr)

        do local_id = 1, grid%nlmax
          trace_epsilon = strain_p((local_id - 1)*SIX_INTEGER + ONE_INTEGER) + &
                          strain_p((local_id - 1)*SIX_INTEGER + TWO_INTEGER) + &
                          strain_p((local_id - 1)*SIX_INTEGER + THREE_INTEGER)
          por_new = por0_p(local_id)/(1.d0 + (1.d0 - por0_p(local_id))*trace_epsilon)
          por_p(local_id) = por_new
        enddo

        call VecRestoreArrayF90(pmc%sim_aux%subsurf_por0, por0_p, ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_strain, strain_p, ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_por, por_p, ierr)

      endif

  end select

end subroutine PMCGeomechanicsSetAuxData

! ************************************************************************** !
!> This routine updates data for geomechanics simulation from other process 
!! models.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 01/01/14
! ************************************************************************** !
subroutine PMCGeomechanicsGetAuxData(this)

  use Option_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"

  class(pmc_geomechanics_type) :: this

  PetscErrorCode :: ierr

  select type(pmc => this)
    class is(pmc_geomechanics_type)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_pres, &
                           pmc%geomech_realization%geomech_field%press, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_pres, &
                         pmc%geomech_realization%geomech_field%press, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%press, &
                            pmc%geomech_realization%geomech_field%press_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%temp, &
                            pmc%geomech_realization%geomech_field%temp_loc, &
                            ONEDOF)

  end select

end subroutine PMCGeomechanicsGetAuxData

end module PMC_Geomechanics_class

#endif
