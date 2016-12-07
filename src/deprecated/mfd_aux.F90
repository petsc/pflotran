module MFD_Aux_module
#include "finclude/petscvec.h"
  use petscvec

  use Connection_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

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
    type(mfd_auxvar_type), pointer :: auxvars(:)
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
!geh: deprecated in PETSc in spring 2014
!    ISLocalToGlobalMapping :: mapping_ltogb_faces ! block form of mapping_ltog
    ISLocalToGlobalMapping :: mapping_ltog_LP  ! petsc vec local to global mapping
!geh: deprecated in PETSc in spring 2014
!    ISLocalToGlobalMapping :: mapping_ltogb_LP ! block form of mapping_ltog
  end type mfd_type

  public :: MFDAuxCreate, MFDAuxDestroy, &
            MFDAuxInit, MFDAuxVarInit, MFDAuxAddFace, & 
            MFDAuxVarDestroy,  &
            MFDAuxInitStiffMatrix, MFDAuxInitResidDerivArrays 

contains

! ************************************************************************** !

function MFDAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  use Option_module

  implicit none
  
  type(mfd_type), pointer :: MFDAuxCreate
  
  type(mfd_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  nullify(aux%auxvars)

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

  MFDAuxCreate => aux
  
end function MFDAuxCreate

! ************************************************************************** !

subroutine MFDAuxInit(aux, num_aux, option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  use Option_module

  implicit none
  
  type(mfd_type) :: aux
  PetscInt :: num_aux, i
  type(option_type) :: option

  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  allocate(aux%auxvars(num_aux))

  do i=1,num_aux
    nullify(aux%auxvars(i)%face_id_gh)
    nullify(aux%auxvars(i)%MassMatrixInv)
    nullify(aux%auxvars(i)%StiffMatrix)
    nullify(aux%auxvars(i)%Rl)
    nullify(aux%auxvars(i)%dRl_dp)
    nullify(aux%auxvars(i)%dRp_dl)
    nullify(aux%auxvars(i)%dRp_dneig)
  end do

end subroutine MFDAuxInit

! ************************************************************************** !

subroutine MFDAuxVarInit(auxvar, numfaces, option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  use Option_module

  implicit none
  
  type(mfd_auxvar_type), pointer :: auxvar
  PetscInt :: numfaces
  type(option_type) :: option

  auxvar%numfaces = numfaces
  if (associated(auxvar%face_id_gh)) deallocate(auxvar%face_id_gh)
  allocate(auxvar%face_id_gh(numfaces))
  auxvar%face_id_gh(1:numfaces) = 0
  
end subroutine MFDAuxVarInit

! ************************************************************************** !

subroutine MFDAuxAddFace(auxvar, option, face_id)
  ! 
  ! MFDAuxVarAddFace: Add face_id to list of faces
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  use Option_module

  implicit none
  
  type(mfd_auxvar_type) :: auxvar
  PetscInt :: face_id
  type(option_type) :: option 

  PetscInt :: i
  PetscBool :: done

  done = PETSC_FALSE

  do i = 1, auxvar%numfaces
     if (auxvar%face_id_gh(i) == 0) then
        auxvar%face_id_gh(i) = face_id
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

subroutine MFDAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  implicit none

  type(mfd_auxvar_type) :: auxvar


  if (associated(auxvar%face_id_gh)) deallocate(auxvar%face_id_gh)
  nullify(auxvar%face_id_gh)

  if (associated(auxvar%MassMatrixInv)) deallocate(auxvar%MassMatrixInv)
  nullify(auxvar%MassMatrixInv)

  if (associated(auxvar%StiffMatrix)) deallocate(auxvar%StiffMatrix)
  nullify(auxvar%StiffMatrix)

  if (associated(auxvar%Rl)) deallocate(auxvar%Rl)
  nullify(auxvar%Rl)

  if (associated(auxvar%dRp_dl)) deallocate(auxvar%dRp_dl)
  nullify(auxvar%dRp_dl)

  if (associated(auxvar%dRl_dp)) deallocate(auxvar%dRl_dp)
  nullify(auxvar%dRl_dp)

  if (associated(auxvar%dRp_dneig)) deallocate(auxvar%dRp_dneig)
  nullify(auxvar%dRp_dneig)

end subroutine MFDAuxVarDestroy

! ************************************************************************** !

subroutine MFDAuxDestroy(aux)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 02/03/10
  ! 

  implicit none

  type(mfd_type), pointer :: aux
  PetscInt :: iaux, sz
  PetscErrorCode :: ierr
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MFDAuxVarDestroy(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  
  if (aux%is_ghosted_local_faces /= 0) then
         call ISDestroy(aux%is_ghosted_local_faces, ierr);CHKERRQ(ierr)
       endif
  if (aux%is_local_local_faces /= 0) then
           call ISDestroy(aux%is_local_local_faces, ierr);CHKERRQ(ierr)
         endif
  if (aux%is_ghosted_petsc_faces /= 0) then
           call ISDestroy(aux%is_ghosted_petsc_faces, ierr);CHKERRQ(ierr)
         endif
  if (aux%is_local_petsc_faces /= 0) then
           call ISDestroy(aux%is_local_petsc_faces, ierr);CHKERRQ(ierr)
         endif
  if (aux%is_ghosts_local_faces /= 0) then
           call ISDestroy(aux%is_ghosts_local_faces, ierr);CHKERRQ(ierr)
         endif
  if (aux%is_ghosts_petsc_faces /= 0) then
           call ISDestroy(aux%is_ghosts_petsc_faces, ierr);CHKERRQ(ierr)
         endif

  if (aux%is_ghosted_local_LP /= 0) then
                 call ISDestroy(aux%is_ghosted_local_LP, ierr);CHKERRQ(ierr)
               endif
  if (aux%is_local_local_LP /= 0)  then
                 call ISDestroy(aux%is_local_local_LP, ierr);CHKERRQ(ierr)
               endif
  if (aux%is_ghosted_petsc_LP /= 0) then 
                 call ISDestroy(aux%is_ghosted_petsc_LP, ierr);CHKERRQ(ierr)
               endif
  if (aux%is_local_petsc_LP /= 0)  then
                 call ISDestroy(aux%is_local_petsc_LP, ierr);CHKERRQ(ierr)
               endif
  if (aux%is_ghosts_local_LP /= 0)  then
                 call ISDestroy(aux%is_ghosts_local_LP, ierr);CHKERRQ(ierr)
               endif

  if (aux%is_ghosts_petsc_LP /= 0)  then
                 call ISDestroy(aux%is_ghosts_petsc_LP, ierr);CHKERRQ(ierr)
               endif

  if (aux%scatter_ltog_faces /= 0) then
                 call VecScatterDestroy(aux%scatter_ltog_faces,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_gtol_faces /= 0) then
                 call VecScatterDestroy(aux%scatter_gtol_faces,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_ltol_faces /= 0) then
                 call VecScatterDestroy(aux%scatter_ltol_faces,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_gtogh_faces /= 0) then
                 call VecScatterDestroy(aux%scatter_gtogh_faces,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_gtol_LP /= 0) then
                 call VecScatterDestroy(aux%scatter_gtol_LP,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_ltog_LP /= 0) then
                 call VecScatterDestroy(aux%scatter_ltog_LP,  &
                                        ierr);CHKERRQ(ierr)
               endif
  if (aux%scatter_ltol_LP /= 0) then
                 call VecScatterDestroy(aux%scatter_ltol_LP,  &
                                        ierr);CHKERRQ(ierr)
               endif

  if (aux%mapping_ltog_faces /= 0) then
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltog_faces,  &
                                                    ierr);CHKERRQ(ierr)
               endif
  
  if (aux%mapping_ltog_LP /= 0) then
                 call ISLocalToGlobalMappingDestroy(aux%mapping_ltog_LP,  &
                                                    ierr);CHKERRQ(ierr)
               endif
    
end subroutine MFDAuxDestroy

! ************************************************************************** !

subroutine MFDAuxInitResidDerivArrays(auxvar, option)

use Option_module

  implicit none

  type(mfd_auxvar_type), pointer :: auxvar
  type(option_type) :: option

  allocate(auxvar%Rl(auxvar%numfaces))
  allocate(auxvar%dRl_dp(auxvar%numfaces))
  allocate(auxvar%dRp_dl(auxvar%numfaces))
  allocate(auxvar%dRp_dneig(auxvar%numfaces))

  auxvar%Rl = 0.
  auxvar%dRl_dp = 0.
  auxvar%dRp_dl = 0.
  auxvar%dRp_dneig = 0.
 
  allocate(auxvar%WB(auxvar%numfaces))
  allocate(auxvar%gr(auxvar%numfaces))

  auxvar%WB = 0.
  auxvar%gr = 0.

end subroutine MFDAuxInitResidDerivArrays

! ************************************************************************** !

subroutine MFDAuxInitStiffMatrix(auxvar, option)

   use Option_module

  implicit none

  type(mfd_auxvar_type), pointer :: auxvar
  type(option_type) :: option

  allocate(auxvar%StiffMatrix(auxvar%numfaces, auxvar%numfaces))

  auxvar%StiffMatrix = 0.

end subroutine MFDAuxInitStiffMatrix


end module MFD_Aux_module
