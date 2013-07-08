module Hydrogeophysics_Factory_module

  use Hydrogeophysics_Simulation_class
  
  implicit none

  private

#include "definitions.h"

  public :: HydrogeophysicsInitialize

contains

! ************************************************************************** !
!
! HydrogeophysicsInitialize: Sets up hydrogeophysics simulation 
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitialize(simulation_base,option)

  use Option_module
  use Hydrogeophysics_Wrapper_module
  use Input_module
  use Simulation_Base_class 
  
  use vars, only : PFE4D_COMM
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: simulation
  PetscMPIInt :: mycolor_mpi, mykey_mpi, pf_e4d_ranks(2)
  PetscInt :: i, num_subsurface_processes, offset
  PetscInt :: local_size
  PetscErrorCode :: ierr

  ! NOTE: PETSc must already have been initialized here!
  if (option%global_commsize < 3) then
    option%io_buffer = 'Must of at least processes allocates to ' // &
      'simulation in order to run hydrogeophysics.'
    call printErrMsg(option)
  endif
  
  num_subsurface_processes = 1
  
  ! split the communicator
  option%mygroup_id = 0
  offset = 0
  if (option%global_rank > num_subsurface_processes-1) then
    option%mygroup_id = 1
    offset = num_subsurface_processes
  endif
  mycolor_mpi = option%mygroup_id
  mykey_mpi = option%global_rank - offset
  call MPI_Comm_split(option%global_comm,mycolor_mpi,mykey_mpi,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)  
  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)
  
  ! create group shared by both master processes
  if (option%myrank == 0) then
    call MPI_Comm_group(option%global_comm,simulation%pf_e4d_grp,ierr)
    pf_e4d_ranks(1) = 0
    pf_e4d_ranks(2) = num_subsurface_processes
    call MPI_Group_incl(option%global_group,TWO_INTEGER_MPI,pf_e4d_ranks, &
                        simulation%pf_e4d_grp,ierr)
    call MPI_Comm_create(option%global_comm,simulation%pf_e4d_grp, &
                         simulation%pf_e4d_comm,ierr)
    call MPI_Comm_rank(simulation%pf_e4d_comm,simulation%pf_e4d_rank,ierr)
    call MPI_Comm_size(simulation%pf_e4d_comm,simulation%pf_e4d_size,ierr)
    PFE4D_COMM = simulation%pf_e4d_comm
  endif  
  
  simulation => HydrogeophysicsCreate(option)
  if (option%mygroup_id == 0) then
    call HydrogeophysicsInitPostPetsc(simulation,option)
    simulation%subsurface_process = PETSC_TRUE
  else
    call HydrogeophysicsWrapperInit(option)
  endif

#if 0  
  ! Petsc vectors are integer quantities, not pointers
  if (simulation%subsurface_process) then
    allocate(int_array(grid%nlmax))
    do local_id = 1, grid%nlmax 
      int_array(local_id) = grid%nG2A(local_id)-1
    enddo
    call ISCreateGeneral(option%mycomm,unstructured_grid%nlmax, &
                         int_array,PETSC_COPY_VALUES,is_tmp,ierr) 
    call VecGetLocalSize(simulation%realization%field%work,local_size,ierr)
  else
    local_size = 0
  endif
  call VecCreate(option%global_comm,simulation%sigma,ierr)
  call VecSetSizes(simulation%sigma,local_size,PETSC_DECIDE,ierr)

  if (simulation%subsurface_process) then
    allocate(int_array(grid%nlmax))
    do local_id = 1, grid%nlmax 
      int_array(local_id) = grid%nG2A(local_id)-1
    enddo
    call ISCreateGeneral(option%global_comm,unstructured_grid%nlmax, &
                         int_array,PETSC_COPY_VALUES,is_global,ierr) 
  else
    call ISCreateGeneral(option%global_comm,ZERO_INTEGER, &
                         int_array,PETSC_COPY_VALUES,is_global,ierr) 
    if (option%myrank == 0) then
      deallocate(int_array)
      call VecGetSize(simulation%sigma,global_size,ierr)
      allocate(int_array(global_size))
      do local_id = 1, global_size
        int_array(local_id) = local_id - 1
      enddo
      call ISCreateGeneral(PETSC_COMM
    endif
  endif
  deallocate(int_array)
  
  
  
  
  deallocate(int_array)
  do local_id = 1, unstructured_grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr)
  call ISDestroy(is_tmp,ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_natural,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  string = 'is_local_natural' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
  !                  PETSC_DETERMINE,vec_tmp,ierr)
  call VecCreate(option%mycomm,vec_tmp,ierr)
  call VecSetSizes(vec_tmp,unstructured_grid%nlmax*ndof,PETSC_DECIDE,ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr)
  call VecSetFromOptions(vec_tmp,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_petsc,vec_tmp, &
                        ugdm%is_local_natural,ugdm%scatter_gton,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_natural,vec_tmp, &
                        ugdm%is_local_petsc,ugdm%scatter_ntog,ierr)
  call VecDestroy(vec_tmp,ierr)  
  
  
  
  
  
  call VecSetFromOptions(simulation%sigma,ierr)
  
#endif  
  
  simulation%hydrogeophysics_coupler%sigma = simulation%sigma
  
  simulation_base => simulation

end subroutine HydrogeophysicsInitialize

! ************************************************************************** !
!
! HydrogeophysicsInitializePostPetsc: Sets up hydrogeophysics simulation 
!                                     framework after to PETSc initialization
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitPostPetsc(simulation, option)

  use Simulation_module
  use Subsurface_Factory_module
  use PMC_Hydrogeophysics_class
  use Option_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  class(hydrogeophysics_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  class(pmc_hydrogeophysics_type), pointer :: hydrogeophysics_coupler
  PetscErrorCode :: ierr
  
  ! Init() is called in SubsurfaceInitializePostPetsc
  call SubsurfaceInitializePostPetsc(simulation, option)
  
  ! add hydrogeophysics coupler to list
  hydrogeophysics_coupler => PMCHydrogeophysicsCreate()
  hydrogeophysics_coupler%option => simulation%option
  hydrogeophysics_coupler%realization => simulation%realization
  simulation%hydrogeophysics_coupler => hydrogeophysics_coupler
  simulation%process_model_coupler_list%below%below => hydrogeophysics_coupler

end subroutine HydrogeophysicsInitPostPetsc

! ************************************************************************** !
!
! HydrogeoInitCommandLineSettings: Initializes hydrogeophysics settings
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeoInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
!  string = '-dummy'
!  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
!  if (option_found) then
!    option%subsurface_simulation_type = dummy
!  endif
  
end subroutine HydrogeoInitCommandLineSettings

end module Hydrogeophysics_Factory_module
