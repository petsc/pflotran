module PMC_Hydrogeophysics_class

  use PMC_Base_class
  use Realization_class
  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(pmc_base_type) :: pmc_hydrogeophysics_type
    class(realization_type), pointer :: realization
    Vec :: solution_seq
    Vec :: solution_mpi
    VecScatter :: pf_to_e4d_scatter
    PetscMPIInt :: pf_to_e4d_master_comm
  contains
    procedure, public :: Init => PMCHydrogeophysicsInit
    procedure, public :: InitializeRun => PMCHydrogeophysicsInitializeRun
    procedure, public :: RunToTime => PMCHydrogeophysicsRunToTime
    procedure, public :: FinalizeRun => PMCHydrogeophysicsFinalizeRun
    procedure, public :: Destroy => PMCHydrogeophysicsDestroy
    procedure, public :: GetAuxData => PMCHydrogeophysicsSynchronize
  end type pmc_hydrogeophysics_type
  
  public :: PMCHydrogeophysicsCreate
  
contains

! ************************************************************************** !

function PMCHydrogeophysicsCreate()
  ! 
  ! Allocates and initializes a new
  ! process_model_coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_hydrogeophysics_type), pointer :: PMCHydrogeophysicsCreate
  
  class(pmc_hydrogeophysics_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCHydrogeophysicsCreate => pmc  
  
end function PMCHydrogeophysicsCreate

! ************************************************************************** !

subroutine PMCHydrogeophysicsInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call PMCBaseInit(this)
  this%name = 'PMCHydrogeophysics'
  nullify(this%realization) 
  this%solution_mpi = 0
  this%solution_seq = 0
  this%pf_to_e4d_scatter = 0

end subroutine PMCHydrogeophysicsInit

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%InitializeRun()')
  
  if (associated(this%below)) then
    call this%below%InitializeRun()
  endif
  
  if (associated(this%next)) then
    call this%next%InitializeRun()
  endif

end subroutine PMCHydrogeophysicsInitializeRun

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Hydrogeophysics_Wrapper_module, only : HydrogeophysicsWrapperStep

  implicit none
  
  class(pmc_hydrogeophysics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  class(pmc_base_type), pointer :: pmc_base
  PetscInt :: local_stop_flag
  
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  call this%GetAuxData()
  
  local_stop_flag = 0
  
  call HydrogeophysicsWrapperStep(sync_time,this%solution_mpi, &
                                  this%solution_seq, &
                                  this%pf_to_e4d_scatter, &
                                  this%pf_to_e4d_master_comm,this%option)

  ! Run neighboring process model couplers
!  if (associated(this%below)) then
!    call this%below%RunToTime(sync_time,local_stop_flag)
!  endif

  ! Run neighboring process model couplers
!  if (associated(this%next)) then
!    call this%next%RunToTime(sync_time,local_stop_flag)
!  endif

  stop_flag = max(stop_flag,local_stop_flag)  
  
end subroutine PMCHydrogeophysicsRunToTime

! ************************************************************************** !

subroutine PMCHydrogeophysicsSynchronize(this)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Variables_module, only : PRIMARY_MOLALITY
  use String_module
!  use Discretization_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

  class(pmc_hydrogeophysics_type) :: this

  PetscReal, pointer :: vec1_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr
  PetscInt, save :: num_calls = 0
  character(len=MAXSTRINGLENGTH) :: filename
  PetscViewer :: viewer
  Vec :: natural_vec
  PetscInt :: i

#if 1
  call RealizationGetVariable(this%realization,this%realization%field%work, &
                              PRIMARY_MOLALITY,ONE_INTEGER,0)
#else
  call DiscretizationCreateVector(this%realization%discretization,ONEDOF, &
                                  natural_vec,NATURAL,this%option)
  if (this%option%myrank == 0) then
    do i = 1, this%realization%patch%grid%nmax
      call VecSetValues(natural_vec,1,i-1,i*1.d0, &
                        INSERT_VALUES,ierr)
    enddo
  endif
  call VecAssemblyBegin(natural_vec,ierr) 
  call VecAssemblyEnd(natural_vec,ierr) 
  call DiscretizationNaturalToGlobal(this%realization%discretization, &
                                      natural_vec, &
                                      this%realization%field%work,ONEDOF)
  call VecDestroy(natural_vec,ierr)
#endif
  call VecGetArrayF90(this%realization%field%work,vec1_ptr,ierr)
  call VecGetArrayF90(this%solution_mpi,vec2_ptr,ierr)
!      vec1_ptr(:) = vec1_ptr(:) + num_calls
  vec2_ptr(:) = vec1_ptr(:)
!  print *, 'PMC update to solution', vec2_ptr(16)
  call VecRestoreArrayF90(this%realization%field%work,vec1_ptr,ierr)
  call VecRestoreArrayF90(this%solution_mpi,vec2_ptr,ierr)

#if 0
  filename = 'pf_solution' // trim(StringFormatInt(num_calls)) // '.txt'
  call PetscViewerASCIIOpen(this%option%mycomm,filename,viewer,ierr)
  call VecView(this%realization%field%work,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  num_calls = num_calls + 1
  
end subroutine PMCHydrogeophysicsSynchronize

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Hydrogeophysics_Wrapper_module, only : HydrogeophysicsWrapperStop
  
  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%FinalizeRun()')
  
  if (this%pf_to_e4d_master_comm /= MPI_COMM_NULL) then
    call HydrogeophysicsWrapperStop(this%option,this%pf_to_e4d_master_comm)
  endif
  
end subroutine PMCHydrogeophysicsFinalizeRun

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsDestroy(this)
  ! 
  ! Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this

  PetscErrorCode :: ierr
  
  call printMsg(this%option,'PMCHydrogeophysics%Destroy()')
  
  nullify(this%realization)
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 

  if (this%solution_seq /= 0) &
    call VecDestroy(this%solution_seq,ierr)
  this%solution_seq = 0
  if (this%pf_to_e4d_scatter /= 0) &
    call VecScatterDestroy(this%pf_to_e4d_scatter, ierr)
  this%pf_to_e4d_scatter = 0
  
end subroutine PMCHydrogeophysicsDestroy
  
end module PMC_Hydrogeophysics_class
