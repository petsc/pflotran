module MFD_module


  use Connection_module
  use Grid_module
  use MFD_Aux_module
  implicit none
  
  private 

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"



  public :: MFDCreateJacobian

contains

! ************************************************************************** !
!
! MFDCreateJacobian: Creates a Jacobian matrix for  the faced based dof
! author: Daniil Svyatskiy 
! date: 05/14/10
!
! ************************************************************************** !


subroutine MFDCreateJacobian(grid, mfd_aux, mat_type, J, option)

 use Option_module
 use Grid_module
 use MFD_Aux_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(grid_type) :: grid
  type(mfd_type) :: mfd_aux
  type(option_type) :: option
  MatType :: mat_type
  Mat :: J


  type(mfd_auxvar_type), pointer :: aux_var
  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt :: icell, iface
  PetscInt :: ghost_face_id, local_face_id, ghost_face_id_n, local_face_id_n
  PetscInt :: icount, jcount
  PetscInt :: ndof_local
  PetscErrorCode :: ierr

  allocate(d_nnz(grid%nlmax_faces))
  allocate(o_nnz(grid%nlmax_faces))

  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0

  do icell = 1, grid%nlmax
    aux_var => MFD_aux%aux_vars(icell)
    do icount = 1, aux_var%numfaces
      ghost_face_id = aux_var%face_id_gh(icount)
      local_face_id = grid%fG2L(ghost_face_id)
      if (local_face_id > 0) then
        do jcount = 1, aux_var%numfaces
          if (icount == jcount) cycle
          ghost_face_id_n = aux_var%face_id_gh(jcount)
          local_face_id_n = grid%fG2L(ghost_face_id_n)
          if (local_face_id_n > 0) then
             d_nnz(local_face_id) = d_nnz(local_face_id) + 1
          else
             o_nnz(local_face_id) = o_nnz(local_face_id) + 1
          end if
        end do
      end if
    end do
  end do

 ! do iface = 1, grid%nlmax_faces
 !   write(*,*) iface, d_nnz(iface), o_nnz(iface)
 ! end do 

  ndof_local = mfd_aux%ndof * grid%nlmax_faces

  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        o_nnz = o_nnz*mfd_aux%ndof
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltogb_faces,ierr)
      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltogb_faces,ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  else
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  endif



  deallocate(d_nnz)
  deallocate(o_nnz)


end subroutine MFDCreateJacobian

end module MFD_module
