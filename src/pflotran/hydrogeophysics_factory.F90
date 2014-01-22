module Hydrogeophysics_Factory_module

  use Hydrogeophysics_Simulation_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

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
  use Input_Aux_module
  use Simulation_Base_class 
  use Discretization_module
  use String_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscviewer.h"
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: simulation
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscInt :: i, num_pflotran_processes, offset, num_slaves
  PetscInt :: local_size
  PetscBool :: option_found
  PetscMPIInt :: mpi_int, process_range(3)
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscInt, allocatable :: int_array(:)
  IS :: is_natural
  IS :: is_petsc
  PetscInt :: istart
  PetscViewer :: viewer
  Vec :: pflotran_solution_vec_mpi, pflotran_solution_vec_seq
  VecScatter :: pflotran_scatter
  
#ifndef E4D
  option%io_buffer = 'Must compile with E4D defined during preprocessing ' // &
    'step in order to use HYDROGEOPHYSICS.'
  call printErrMsg(option)
#endif

  string = '-num_slaves'
  num_slaves = -999
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) num_slaves = i

  ! NOTE: PETSc must already have been initialized here!
  if (option%global_commsize < 3) then
    option%io_buffer = 'At least 3 processes must be allocated to ' // &
      'simulation in order to run hydrogeophysics.'
    call printErrMsg(option)
  endif

  simulation => HydrogeophysicsCreate(option)
  
  if (num_slaves < -998) then
    num_pflotran_processes = option%global_commsize / 2
    num_slaves = option%global_commsize - num_pflotran_processes - 1
  else if (num_slaves <= 0) then
    option%io_buffer = 'Number of slaves must be greater than zero. ' // &
      'Currently set to ' // StringFormatInt(num_slaves) // '.'
    call printErrMsg(option)
  else
    if (num_slaves > option%global_commsize - 2) then
      option%io_buffer = 'Too many slave processes allocated to ' // &
        'simulation: ' // StringFormatInt(num_slaves)
      call printErrMsg(option)
    endif
    num_pflotran_processes = option%global_commsize - num_slaves - 1
  endif
  
  write(option%io_buffer,*) 'Number of E4D processes: ', &
    StringFormatInt(num_slaves+1)
  call printMsg(option)
  write(option%io_buffer,*) 'Number of PFLOTRAN processes: ', &
    StringFormatInt(num_pflotran_processes)
  call printMsg(option)
  
  ! split the communicator
  option%mygroup_id = 0
  offset = 0
  if (option%global_rank > num_pflotran_processes-1) then
    option%mygroup_id = 1
    offset = num_pflotran_processes
  endif

  if (option%mygroup_id == 0) then
    simulation%pflotran_process = PETSC_TRUE
  endif

  mycolor_mpi = option%mygroup_id
  mykey_mpi = option%global_rank - offset
  call MPI_Comm_split(option%global_comm,mycolor_mpi,mykey_mpi,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)  
  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)
  
  ! create group shared by both master processes
  call MPI_Comm_group(option%global_comm,option%global_group,ierr)
  call MPI_Group_size(option%mygroup,mpi_int,ierr)
!print *, 1, option%global_rank, option%global_group, option%global_comm
!print *, 2, option%global_rank, option%myrank, option%mygroup, option%mycomm, mpi_int
  mpi_int = 1
  process_range(1) = 0
  process_range(2) = num_pflotran_processes ! includes e4d master due to 
  process_range(3) = 1                        ! zero-based indexing
  call MPI_Group_range_incl(option%global_group,mpi_int,process_range, &
                            simulation%pf_e4d_scatter_grp,ierr)
!print *, 3, option%global_rank, simulation%pf_e4d_scatter_grp
  call MPI_Comm_create(option%global_comm,simulation%pf_e4d_scatter_grp, &
                       simulation%pf_e4d_scatter_comm,ierr)
!print *, 4, option%global_rank, simulation%pf_e4d_scatter_comm, simulation%pf_e4d_master_comm
  if (simulation%pf_e4d_scatter_comm /= MPI_COMM_NULL) then
    call MPI_Comm_rank(simulation%pf_e4d_scatter_comm, &
                       simulation%pf_e4d_scatter_rank,ierr)
    call MPI_Comm_size(simulation%pf_e4d_scatter_comm, &
                       simulation%pf_e4d_scatter_size,ierr)
!    PFE4D_COMM = simulation%pf_e4d_scatter_comm
!print *, 5, option%global_rank, simulation%pf_e4d_scatter_rank, simulation%pf_e4d_scatter_size
    ! remove processes between pf_master and e4d_master
    mpi_int = 1
    process_range(1) = 1
    process_range(2) = num_pflotran_processes-1
    process_range(3) = 1
    ! if there are no process ranks to remove, set mpi_int to zero
    if (process_range(2) - process_range(1) < 0) mpi_int = 0
    call MPI_Group_range_excl(simulation%pf_e4d_scatter_grp,mpi_int, &
                              process_range,simulation%pf_e4d_master_grp,ierr)
!print *, 6, option%global_rank, simulation%pf_e4d_master_grp, ierr
    call MPI_Comm_create(simulation%pf_e4d_scatter_comm, &
                         simulation%pf_e4d_master_grp, &
                         simulation%pf_e4d_master_comm,ierr)
!print *, 7, option%global_rank, simulation%pf_e4d_master_comm
    if (simulation%pf_e4d_master_comm /= MPI_COMM_NULL) then
      call MPI_Comm_rank(simulation%pf_e4d_master_comm, &
                         simulation%pf_e4d_master_rank,ierr)
      call MPI_Comm_size(simulation%pf_e4d_master_comm, &
                         simulation%pf_e4d_master_size,ierr)
!      PFE4D_MASTER_COMM = simulation%pf_e4d_master_comm
!print *, 8, option%global_rank, simulation%pf_e4d_master_rank, simulation%pf_e4d_master_size
    endif
  endif

!print *, 9, option%global_rank, simulation%pf_e4d_scatter_comm, simulation%pf_e4d_master_comm, MPI_COMM_NULL
  
  if (simulation%pflotran_process) then
    call HydrogeophysicsInitPostPetsc(simulation,option)
  else
    option%io_rank = -1 ! turn off I/O from E4D processes.
  endif

!#define DEBUG

  !   PFLOTRAN subsurface processes      E4D master process
  if (simulation%pflotran_process .or. option%myrank == 0) then 
    if (simulation%pflotran_process) then
      simulation%hydrogeophysics_coupler%pf_to_e4d_master_comm = &
        simulation%pf_e4d_master_comm
    endif

    ! create mpi Vec that includes all PFLOTRAN processes and the E4D master
    call VecCreate(simulation%pf_e4d_scatter_comm,pflotran_solution_vec_mpi, &
                   ierr)
    if (simulation%pflotran_process) then
      call VecGetLocalSize(simulation%realization%field%work,local_size,ierr)
    else ! E4D master process
      local_size = 0
    endif
    call VecSetSizes(pflotran_solution_vec_mpi,local_size,PETSC_DECIDE,ierr)
    call VecSetFromOptions(pflotran_solution_vec_mpi,ierr)

    allocate(int_array(local_size))
    if (simulation%pflotran_process) then
      int_array = 0
      do i = 1, local_size
        int_array(i) =  simulation%realization%patch%grid%nG2A( &
                          simulation%realization%patch%grid%nL2G(i))

      enddo
print *, option%myrank, int_array
      int_array = int_array - 1
    endif
    call ISCreateGeneral(simulation%pf_e4d_scatter_comm,local_size,int_array, &
!    call ISCreateGeneral(PETSC_COMM_SELF,local_size,int_array, &
                         PETSC_COPY_VALUES,is_natural,ierr)
    deallocate(int_array)
    if (simulation%pflotran_process) then
      call VecGetOwnershipRange(simulation%realization%field%work,istart, &
                                PETSC_NULL_INTEGER,ierr)
    endif
!    istart = 0
    call ISCreateStride(simulation%pf_e4d_scatter_comm,local_size,istart,1, &
                        is_petsc,ierr)

#ifdef DEBUG
  string = 'is_petsc.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer,ierr)
  call ISView(is_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  string = 'is_natural.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer,ierr)
!  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call ISView(is_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif


    ! create seq Vec on each (all PFLOTRAN processes and the E4D master)
    call VecGetSize(pflotran_solution_vec_mpi,local_size,ierr)
    if (simulation%pflotran_process) local_size = 0
    ! make global size the local size of pflotran_solution_vec_seq on E4D master
!    call VecCreate(PETSC_COMM_SELF,pflotran_solution_vec_seq,ierr)
    call VecCreate(simulation%pf_e4d_scatter_comm,pflotran_solution_vec_seq,ierr)
    call VecSetSizes(pflotran_solution_vec_seq,local_size,PETSC_DECIDE,ierr)
    call VecSetFromOptions(pflotran_solution_vec_seq,ierr)

#ifdef DEBUG
  string = 'pflotran_solution_vec_mpi.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer,ierr)
  call VecView(pflotran_solution_vec_mpi,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  string = 'pflotran_solution_vec_seq' // trim(StringFormatInt(option%global_rank)) // '.txt'
!  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer,ierr)
  call VecView(pflotran_solution_vec_seq,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
    ! create scatter between mpi Vec and local seq Vecs (only E4D master is 
    ! relevant)
    call VecScatterCreate(pflotran_solution_vec_mpi,is_petsc, &
                          pflotran_solution_vec_seq,is_natural, &
                          pflotran_scatter,ierr)
#ifdef DEBUG
  string = 'scatter.txt'
  call PetscViewerASCIIOpen(simulation%pf_e4d_scatter_comm,trim(string),viewer,ierr)
  call VecScatterView(pflotran_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
    call ISDestroy(is_natural,ierr)
    call ISDestroy(is_petsc,ierr)
  endif
!print *, 'End  -----------'

  simulation_base => simulation

  if (simulation%pflotran_process) then
    simulation%hydrogeophysics_coupler%solution_mpi = pflotran_solution_vec_mpi
    simulation%hydrogeophysics_coupler%solution_seq = pflotran_solution_vec_seq
    simulation%solution_mpi = pflotran_solution_vec_mpi
    simulation%hydrogeophysics_coupler%pf_to_e4d_scatter = pflotran_scatter
    simulation%hydrogeophysics_coupler%pf_to_e4d_master_comm = &
      simulation%pf_e4d_master_comm
  else
    call HydrogeophysicsWrapperInit(option, &
                                    pflotran_solution_vec_mpi, &
                                    pflotran_solution_vec_seq, &
                                    pflotran_scatter, &
                                    simulation%pf_e4d_master_comm)
  endif
#undef DEBUG

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
  if (associated(simulation%process_model_coupler_list%below)) then
    simulation%process_model_coupler_list%below%below => hydrogeophysics_coupler
  else
    simulation%process_model_coupler_list%below => hydrogeophysics_coupler
  endif

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
  use Input_Aux_module
  
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
