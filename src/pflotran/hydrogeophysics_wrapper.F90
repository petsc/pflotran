module Hydrogeophysics_Wrapper_module
 
  implicit none

  private

#include "definitions.h"

  public :: HydrogeophysicsWrapperInit, &
            HydrogeophysicsWrapperStart, &
            HydrogeophysicsWrapperStep, &
            HydrogeophysicsWrapperStop, &
            HydrogeophysicsWrapperDestroy

contains

! ************************************************************************** !
!
! HydrogeophysicsWrapperInit: Initializes the hydrogeophysics module
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperInit(option, &
                                      pflotran_solution_vec_mpi_, &
                                      pflotran_solution_vec_seq_, &
                                      pflotran_scatter_, &
                                      pf_e4d_master_comm)
  
  use Option_module

#ifdef E4D
  use vars, only : E4D_COMM, my_rank, n_rank, PFE4D_MASTER_COMM, &
                   pflotran_solution_vec_mpi, pflotran_solution_vec_seq, &
                   pflotran_scatter
  use e4d_setup, only : setup_e4d
  use e4d_run, only: run_e4D
#endif
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type) :: option
  Vec :: pflotran_solution_vec_mpi_
  Vec :: pflotran_solution_vec_seq_
  VecScatter :: pflotran_scatter_
  PetscMPIInt :: pf_e4d_master_comm
  
!  call printMsg(option,'HydrogeophysicsWrapperInit()')
  
  ! from e4d_vars.F90
#ifdef E4D
  E4D_COMM = option%mycomm
  PFE4D_MASTER_COMM = pf_e4d_master_comm
  my_rank = option%myrank
  n_rank = option%mycommsize
  pflotran_solution_vec_mpi = pflotran_solution_vec_mpi_
  pflotran_solution_vec_seq = pflotran_solution_vec_seq_
  pflotran_scatter = pflotran_scatter_

  call setup_e4d
  call run_e4d
#endif
  
end subroutine HydrogeophysicsWrapperInit

! ************************************************************************** !
!
! HydrogeophysicsWrapperStart: Starts the hydrogeophysics forward simulation 
!                              loop
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStart(option)
  
  use Option_module

  implicit none

  type(option_type) :: option
  
!  call printMsg(option,'HydrogeophysicsWrapperStart()')

end subroutine HydrogeophysicsWrapperStart

! ************************************************************************** !
!
! HydrogeophysicsWrapperStep: Performs a forward simulation
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStep(solution_mpi,solution_seq,scatter, &
                                      comm,option)
  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec :: solution_mpi
  Vec :: solution_seq
  VecScatter :: scatter
  PetscMPIInt :: comm
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
!  call printMsg(option,'HydrogeophysicsWrapperStep()')
  
  ! Bcasting a 1 to E4D tells it to run a forward simulation 
  if (comm /= MPI_COMM_NULL) then
    call MPI_Bcast(ONE_INTEGER_MPI,ONE_INTEGER_MPI,MPI_INTEGER, &
                   ZERO_INTEGER_MPI,comm,ierr)
  endif
  call VecScatterBegin(scatter,solution_mpi,solution_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,solution_mpi,solution_seq, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine HydrogeophysicsWrapperStep

! ************************************************************************** !
!
! HydrogeophysicsWrapperStop: Stops the hydrogeophysics forward simulation 
!                             loop
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStop(option,comm)
  
  use Option_module

  implicit none

  type(option_type) :: option
  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
!  call printMsg(option,'HydrogeophysicsWrapperStop()')
  
  ! Bcasting a 0 to E4D tells it to stop
  call MPI_Bcast(ZERO_INTEGER_MPI,ONE_INTEGER_MPI,MPI_INTEGER, &
                 ZERO_INTEGER_MPI,comm,ierr)

end subroutine HydrogeophysicsWrapperStop

! ************************************************************************** !
!
! HydrogeophysicsWrapperDestroy: Destroys the contents of the hydrogeophysics 
!                                module
! author: Glenn Hammond
! date: 07/02/13
!
!
! ************************************************************************** !
recursive subroutine HydrogeophysicsWrapperDestroy(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  
!  call printMsg(option,'HydrogeophysicsWrapperDestroy()')

end subroutine HydrogeophysicsWrapperDestroy

end module Hydrogeophysics_Wrapper_module
