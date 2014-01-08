module MFD_Aux_module


  use Connection_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"


  type, public :: mfd_auxvar_type
    PetscInt :: numfaces
    PetscInt, pointer :: face_id_gh(:)
    PetscReal, pointer :: MassMatrixInv(:,:)
    PetscReal, pointer :: StiffMatrix(:,:)
    PetscReal :: Rp, dRp_dp
    PetscReal, pointer :: Rl(:), dRp_dl(:), dRl_dp(:), dRp_dneig(:)
    PetscReal, pointer :: WB(:), gr(:)  


  end type mfd_auxvar_type
  
  type, public :: mfd_type
    PetscInt :: num_aux
    PetscInt :: ndof
    type(mfd_auxvar_type), pointer :: aux_vars(:)
    IS :: is_ghosted_local_faces ! IS for ghosted faces with local on-processor numbering
    IS :: is_local_local_faces ! IS for local faces with local on-processor numbering
    IS :: is_ghosted_petsc_faces ! IS for ghosted faces with petsc numbering
    IS :: is_local_petsc_faces ! IS for local faces with petsc numbering
    IS :: is_ghosts_local_faces ! IS for ghosted faces with local on-processor numbering
    IS :: is_ghosts_petsc_faces ! IS for ghosted faces with petsc numbering
    VecScatter :: scatter_ltog_faces ! scatter context for local to global updates
    VecScatter :: scatter_gtol_faces ! scatter context for global to local updates
    VecScatter :: scatter_ltol_faces ! scatter context for local to local updates
    VecScatter :: scatter_gtogh_faces ! scatter context for local to local updates
    IS :: is_ghosted_local_LP ! IS for ghosted faces with local on-processor numbering
    IS :: is_local_local_LP ! IS for local faces with local on-processor numbering
    IS :: is_ghosted_petsc_LP ! IS for ghosted faces with petsc numbering
    IS :: is_local_petsc_LP ! IS for local faces with petsc numbering
    IS :: is_ghosts_local_LP ! IS for ghosted faces with local on-processor numbering
    IS :: is_ghosts_petsc_LP ! IS for ghosted faces with petsc numbering
    VecScatter :: scatter_gtol_LP ! scatter context for global to local updates
    VecScatter :: scatter_ltog_LP ! scatter context for global to local updates
    VecScatter :: scatter_ltol_LP ! scatter context for local to local updates
    ISLocalToGlobalMapping :: mapping_ltog_faces  ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb_faces ! block form of mapping_ltog
    ISLocalToGlobalMapping :: mapping_ltog_LP  ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb_LP ! block form of mapping_ltog
  end type mfd_type

  public :: MFDAuxCreate, MFDAuxDestroy, &
            MFDAuxInit, MFDAuxVarInit, MFDAuxAddFace, & 
            MFDAuxVarDestroy,  &
            MFDAuxInitStiffMatrix, MFDAuxInitResidDerivArrays 

contains


! ************************************************************************** !
!
! MFDAuxCreate: Allocate and initialize auxiliary object
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

  aux%is_ghosted_local_LP  = 0
  aux%is_local_local_LP  = 0
  aux%is_ghosted_petsc_LP = 0
  aux%is_local_petsc_LP = 0
  aux%is_ghosts_local_LP  = 0
  aux%is_ghosts_petsc_LP  = 0
        
  aux%scatter_gtol_LP = 0
  aux%scatter_ltog_LP = 0
  aux%scatter_ltol_LP = 0

  aux%mapping_ltog_LP  = 0  
  aux%mapping_ltogb_LP = 0 

  MFDAuxCreate => aux
  
end function MFDAuxCreate

! ************************************************************************** !
!
! MFDAuxInit: Initialize auxiliary object
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
    nullify(aux%aux_vars(i)%Rl)
    nullify(aux%aux_vars(i)%dRl_dp)
    nullify(aux%aux_vars(i)%dRp_dl)
    nullify(aux%aux_vars(i)%dRp_dneig)
  end do

end subroutine MFDAuxInit

! ************************************************************************** !
!
! MFDAuxVarInit: Initialize auxiliary object
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
  PetscBool :: done

  done = PETSC_FALSE

  do i = 1, aux_var%numfaces
     if (aux_var%face_id_gh(i) == 0) then
        aux_var%face_id_gh(i) = face_id
        done = PETSC_TRUE
        exit
     endif
  enddo

  if (.not.done) then
     call printMsg(option, "Imposible to add face to  MFDAuxVar")
     stop
  endif
  
end subroutine MFDAuxAddFace
  
! ************************************************************************** !
!
! MFDAuxVarDestroy: Deallocates a mode auxiliary object
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

  if (associated(aux_var%Rl)) deallocate(aux_var%Rl)
  nullify(aux_var%Rl)

  if (associated(aux_var%dRp_dl)) deallocate(aux_var%dRp_dl)
  nullify(aux_var%dRp_dl)

  if (associated(aux_var%dRl_dp)) deallocate(aux_var%dRl_dp)
  nullify(aux_var%dRl_dp)

  if (associated(aux_var%dRp_dneig)) deallocate(aux_var%dRp_dneig)
  nullify(aux_var%dRp_dneig)

end subroutine MFDAuxVarDestroy

! ************************************************************************** !
!
! MFDAuxDestroy: Deallocates a mode auxiliary object
! author: Daniil Svyatskiy
! date: 02/03/10
!
! ************************************************************************** !
subroutine MFDAuxDestroy(aux)

  implicit none

  type(mfd_type), pointer :: aux
  PetscInt :: iaux, sz
  PetscErrorCode :: ierr
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call MFDAuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
  
  if (aux%is_ghosted_local_faces /= 0) &
         call ISDestroy(aux%is_ghosted_local_faces, ierr)
  if (aux%is_local_local_faces /= 0) &
           call ISDestroy(aux%is_local_local_faces, ierr)
  if (aux%is_ghosted_petsc_faces /= 0) &
           call ISDestroy(aux%is_ghosted_petsc_faces, ierr)
  if (aux%is_local_petsc_faces /= 0) &
           call ISDestroy(aux%is_local_petsc_faces, ierr)
  if (aux%is_ghosts_local_faces /= 0) &
           call ISDestroy(aux%is_ghosts_local_faces, ierr)
  if (aux%is_ghosts_petsc_faces /= 0) &
           call ISDestroy(aux%is_ghosts_petsc_faces, ierr)

  if (aux%is_ghosted_local_LP /= 0) &
                 call ISDestroy(aux%is_ghosted_local_LP, ierr)
  if (aux%is_local_local_LP /= 0)  &
                 call ISDestroy(aux%is_local_local_LP, ierr)
  if (aux%is_ghosted_petsc_LP /= 0) & 
                 call ISDestroy(aux%is_ghosted_petsc_LP, ierr)
  if (aux%is_local_petsc_LP /= 0)  &
                 call ISDestroy(aux%is_local_petsc_LP, ierr) 
  if (aux%is_ghosts_local_LP /= 0)  &
                 call ISDestroy(aux%is_ghosts_local_LP, ierr)

  if (aux%is_ghosts_petsc_LP /= 0)  &
                 call ISDestroy(aux%is_ghosts_petsc_LP, ierr)

  if (aux%scatter_ltog_faces /= 0) &
                 call VecScatterDestroy(aux%scatter_ltog_faces, ierr)
  if (aux%scatter_gtol_faces /= 0) &
                 call VecScatterDestroy(aux%scatter_gtol_faces, ierr)
  if (aux%scatter_ltol_faces /= 0) &
                 call VecScatterDestroy(aux%scatter_ltol_faces, ierr)
  if (aux%scatter_gtogh_faces /= 0) &
                 call VecScatterDestroy(aux%scatter_gtogh_faces, ierr)
  if (aux%scatter_gtol_LP /= 0) &
                 call VecScatterDestroy(aux%scatter_gtol_LP, ierr)
  if (aux%scatter_ltog_LP /= 0) &
                 call VecScatterDestroy(aux%scatter_ltog_LP, ierr)
  if (aux%scatter_ltol_LP /= 0) &
                 call VecScatterDestroy(aux%scatter_ltol_LP, ierr)

  if (aux%mapping_ltog_faces /= 0) &
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltog_faces, ierr) 
  
  if (aux%mapping_ltogb_faces /= 0) &
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltogb_faces , ierr) 
  
  if (aux%mapping_ltog_LP /= 0) &
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltog_LP, ierr) 
  
  if (aux%mapping_ltogb_LP /= 0) &
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltogb_LP , ierr) 
  
    
end subroutine MFDAuxDestroy

! ************************************************************************** !
! ************************************************************************** !
subroutine MFDAuxInitResidDerivArrays(aux_var, option)

use Option_module

  implicit none

  type(mfd_auxvar_type), pointer :: aux_var
  type(option_type) :: option

  allocate(aux_var%Rl(aux_var%numfaces))
  allocate(aux_var%dRl_dp(aux_var%numfaces))
  allocate(aux_var%dRp_dl(aux_var%numfaces))
  allocate(aux_var%dRp_dneig(aux_var%numfaces))

  aux_var%Rl = 0.
  aux_var%dRl_dp = 0.
  aux_var%dRp_dl = 0.
  aux_var%dRp_dneig = 0.
 
  allocate(aux_var%WB(aux_var%numfaces))
  allocate(aux_var%gr(aux_var%numfaces))

  aux_var%WB = 0.
  aux_var%gr = 0.

end subroutine MFDAuxInitResidDerivArrays

! ************************************************************************** !
! ************************************************************************** !
subroutine MFDAuxInitStiffMatrix(aux_var, option)

   use Option_module

  implicit none

  type(mfd_auxvar_type), pointer :: aux_var
  type(option_type) :: option

  allocate(aux_var%StiffMatrix(aux_var%numfaces, aux_var%numfaces))

  aux_var%StiffMatrix = 0.

end subroutine MFDAuxInitStiffMatrix


end module MFD_Aux_module
