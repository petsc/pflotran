module MFD_Aux_module


  use Connection_module
  implicit none
  
  private 

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"


  type, public :: mfd_auxvar_type
    PetscInt :: numfaces
    PetscInt, pointer :: face_id_gh(:)
    PetscScalar, pointer :: MassMatrixInv(:,:)
    PetscScalar, pointer :: StiffMatrix(:,:)
  end type mfd_auxvar_type
  
  type, public :: mfd_type
    PetscInt :: num_aux
    PetscInt :: ndof
    type(mfd_auxvar_type), pointer :: aux_vars(:)
    IS :: is_ghosted_local_faces ! IS for ghosted cells with local on-processor numbering
    IS :: is_local_local_faces ! IS for local cells with local on-processor numbering
    IS :: is_ghosted_petsc_faces ! IS for ghosted cells with petsc numbering
    IS :: is_local_petsc_faces ! IS for local cells with petsc numbering
    IS :: is_ghosts_local_faces ! IS for ghosted cells with local on-processor numbering
    IS :: is_ghosts_petsc_faces ! IS for ghosted cells with petsc numbering
    VecScatter :: scatter_ltog_faces ! scatter context for local to global updates
    VecScatter :: scatter_gtol_faces ! scatter context for global to local updates
    VecScatter :: scatter_ltol_faces ! scatter context for local to local updates
    VecScatter :: scatter_gtogh_faces ! scatter context for local to local updates
    ISLocalToGlobalMapping :: mapping_ltog_faces  ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb_faces ! block form of mapping_ltog
  end type mfd_type

  public :: MFDAuxCreate, MFDAuxDestroy, &
            MFDAuxInit, MFDAuxVarInit, MFDAuxAddFace, & 
            MFDAuxVarDestroy, MFDAuxGenerateMassMatrixInv 

contains


! ************************************************************************** !
!
! MFDAuxCreate: Allocate and initialize auxilliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
function MFDAuxCreate()

  use Option_module

  implicit none
  
  type(mfd_type), pointer :: MFDAuxCreate
  
  type(mfd_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  nullify(aux%aux_vars)

  aux%is_ghosted_local_faces = 0
 
  aux%is_local_local_faces = 0 
  aux%is_ghosted_petsc_faces = 0
  aux%is_local_petsc_faces  = 0
  aux%is_ghosts_local_faces = 0
  aux%is_ghosts_petsc_faces = 0
  aux%scatter_ltog_faces = 0
  aux%scatter_gtol_faces = 0
  aux%scatter_ltol_faces = 0
  aux%scatter_gtogh_faces = 0
  aux%mapping_ltog_faces = 0
  aux%mapping_ltogb_faces = 0

 
  MFDAuxCreate => aux
  
end function MFDAuxCreate





! ************************************************************************** !
!
! MFDAuxInit: Initialize auxilliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxInit(aux, num_aux, option)

  use Option_module

  implicit none
  
  type(mfd_type) :: aux
  PetscInt :: num_aux, i
  type(option_type) :: option

  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  allocate(aux%aux_vars(num_aux))

  do i=1,num_aux
    nullify(aux%aux_vars(i)%face_id_gh)
    nullify(aux%aux_vars(i)%MassMatrixInv)
    nullify(aux%aux_vars(i)%StiffMatrix)
  end do


  
end subroutine MFDAuxInit


! ************************************************************************** !
!
! MFDAuxVarInit: Initialize auxilliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxVarInit(aux_var, numfaces, option)

  use Option_module

  implicit none
  
  type(mfd_auxvar_type), pointer :: aux_var
  PetscInt :: numfaces
  type(option_type) :: option

  aux_var%numfaces = numfaces
  if (associated(aux_var%face_id_gh)) deallocate(aux_var%face_id_gh)
  allocate(aux_var%face_id_gh(numfaces))
  aux_var%face_id_gh(1:numfaces) = 0


  
end subroutine MFDAuxVarInit



! ************************************************************************** !
!
! MFDAuxVarAddFace: Add face_id to list of faces
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxAddFace(aux_var, option, face_id)

  use Option_module

  implicit none
  
  type(mfd_auxvar_type) :: aux_var
  PetscInt :: face_id
  type(option_type) :: option 

  PetscInt :: i
  logical :: done

  done = 0

  do i = 1, aux_var%numfaces
     if (aux_var%face_id_gh(i) == 0) then
        aux_var%face_id_gh(i) = face_id
        done = 1
        exit
     endif
  enddo

  if (done.eqv..FALSE.) then
     call printMsg(option, "Imposible to add face to  MFDAuxVar")
     stop
  endif

  
end subroutine MFDAuxAddFace

  
! ************************************************************************** !
!
! MFDAuxVarDestroy: Deallocates a mode auxilliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxVarDestroy(aux_var)

  implicit none

  type(mfd_auxvar_type) :: aux_var


  if (associated(aux_var%face_id_gh)) deallocate(aux_var%face_id_gh)
  nullify(aux_var%face_id_gh)

  if (associated(aux_var%MassMatrixInv)) deallocate(aux_var%MassMatrixInv)
  nullify(aux_var%MassMatrixInv)

  if (associated(aux_var%StiffMatrix)) deallocate(aux_var%StiffMatrix)
  nullify(aux_var%StiffMatrix)

  

end subroutine MFDAuxVarDestroy

! ************************************************************************** !
!
! MFDAuxDestroy: Deallocates a mode auxilliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxDestroy(aux)

  implicit none

  type(mfd_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call MFDAuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
    
end subroutine MFDAuxDestroy

! ************************************************************************** !
!
! MFDAuxGenerateMassMatrixInv: Create a mass matrix for cell
! author: Daniil Svyatskiy
! date: 05/25/10
!
! ************************************************************************** !
subroutine MFDAuxGenerateMassMatrixInv(aux_var, volume, PermTensor, option)

 use Option_module

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

  type(mfd_auxvar_type), pointer :: aux_var
  PetscScalar :: PermTensor(3,3)
  PetscScalar :: volume
  type(option_type) :: option


  PetscInt :: iface

  if (volume==0.) then
    option%io_buffer = 'Cell volume iz zero'
    call printErrMsg(option)
  end if


  allocate(aux_var%MassMatrixInv(aux_var%numfaces, aux_var%numfaces))
  allocate(aux_var%StiffMatrix(aux_var%numfaces, aux_var%numfaces))

  aux_var%MassMatrixInv = 0.

  do iface = 1, aux_var%numfaces
    aux_var%MassMatrixInv(iface, iface) = (2./volume)*PermTensor(1,1)
  end do 

end subroutine MFDAuxGenerateMassMatrixInv



end module MFD_Aux_module
