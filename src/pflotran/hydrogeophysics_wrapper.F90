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
subroutine HydrogeophysicsWrapperInit(option)
  
  use Option_module

  use vars
  
  implicit none

  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperInit()')
  
  ! from vars
!  integer, dimension(:), allocatable :: e4d_ranks,pf_e4d_ranks
!  integer :: mpi_comm_grp,mpi_e4d_grp,mpi_pfe4d_grp,i
!  integer :: my_wrank,my_pfe4d_rank,n_pfe4drank             !!my mpi rank
!  integer :: tn_rank  
  
  my_wrank = -999
  tn_rank = -999
  
  mpi_comm_grp = -999
  mpi_e4d_grp = -999
  mpi_pfe4d_grp = -999
  my_pfe4d_rank = -999
  n_pfe4drank = -999
  
  E4D_COMM = option%mycomm
  my_rank = option%myrank
  n_rank = option%mycommsize
  
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
  
  call printMsg(option,'HydrogeophysicsWrapperStart()')
  
end subroutine HydrogeophysicsWrapperStart

! ************************************************************************** !
!
! HydrogeophysicsWrapperStep: Performs a forward simulation
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStep(sigma,option)
  
  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec :: sigma
  type(option_type) :: option
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call printMsg(option,'HydrogeophysicsWrapperStep()')
  
  call VecGetArrayF90(sigma,vec_ptr,ierr)
  write(option%io_buffer,'(es13.6)') vec_ptr(365)
  call printMsg(option)
  call VecRestoreArrayF90(sigma,vec_ptr,ierr)
  
#if 0  
  call MPI_Bcast(ONE_INTEGER_MPI,ONE_INTEGER_MPI,MPI_INTEGER, &
                 ZERO_INTEGER_MPI,simulation%pf_e4d_comm,ierr)
  call VecScatter()
#endif  
  
end subroutine HydrogeophysicsWrapperStep

! ************************************************************************** !
!
! HydrogeophysicsWrapperStop: Stops the hydrogeophysics forward simulation 
!                             loop
! author: Glenn Hammond
! date: 07/02/13
!
! ************************************************************************** !
subroutine HydrogeophysicsWrapperStop(option)
  
  use Option_module

  implicit none

  type(option_type) :: option
  
  call printMsg(option,'HydrogeophysicsWrapperStop()')
  
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
  
  call printMsg(option,'HydrogeophysicsWrapperDestroy()')

end subroutine HydrogeophysicsWrapperDestroy

end module Hydrogeophysics_Wrapper_module
