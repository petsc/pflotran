module Matrix_Block_Aux_module

  ! this module cannot depend on any other modules beside Option_module

  implicit none
  
  private 

#include "definitions.h"
 
  type, public :: matrix_block_auxvar_type

    PetscReal, pointer :: dtotal(:,:,:)
    
  end type matrix_block_auxvar_type
  
  type, public :: matrix_block_info_type
    PetscInt :: dim1
    PetscInt :: dim2
    PetscInt :: dim3
  end type matrix_block_info_type
  
  interface MatrixBlockAuxVarInit
    module procedure MatrixBlockAuxVarInit1
    module procedure MatrixBlockAuxVarInit2
  end interface MatrixBlockAuxVarInit
  
  public :: MatrixBlockAuxVarCreate, &
            MatrixBlockAuxVarInit, &
            MatrixBlockAuxVarCopy, &
            MatrixBlockAuxVarDestroy, &
            MatrixBlockInfoCreate, &
            MatrixBlockInfoDestroy
            
contains

! ************************************************************************** !
!
! MatrixBlockAuxCreate: Allocate and initialize auxiliary object
! author: Glenn Hammond
! date: 03/04/2010
!
! ************************************************************************** !
function MatrixBlockAuxVarCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(matrix_block_auxvar_type), pointer :: MatrixBlockAuxVarCreate
  
  type(matrix_block_auxvar_type), pointer :: aux

  allocate(aux)  
  nullify(aux%dtotal)

  MatrixBlockAuxVarCreate => aux
  
end function MatrixBlockAuxVarCreate

! ************************************************************************** !
!
! MatrixBlockAuxVarInit1: Initialize auxiliary object
! author: Glenn Hammond
! date: 03/04/2010
!
! ************************************************************************** !
subroutine MatrixBlockAuxVarInit1(aux_var,dim1,dim2,dim3,option)

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: aux_var
  type(matrix_block_info_type) :: matrix_info
  PetscInt :: dim1
  PetscInt :: dim2
  PetscInt :: dim3
  type(option_type) :: option  
  
  allocate(aux_var%dtotal(dim1,dim2,dim3))
  aux_var%dtotal = 0.d0
  
end subroutine MatrixBlockAuxVarInit1

! ************************************************************************** !
!
! MatrixBlockAuxVarInit2: Initialize auxiliary object
! author: Glenn Hammond
! date: 03/04/2010
!
! ************************************************************************** !
subroutine MatrixBlockAuxVarInit2(aux_var,matrix_info,option)

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: aux_var
  type(matrix_block_info_type) :: matrix_info
  type(option_type) :: option  
  
  allocate(aux_var%dtotal(matrix_info%dim1,matrix_info%dim2,matrix_info%dim3))
  aux_var%dtotal = 0.d0
  
end subroutine MatrixBlockAuxVarInit2

! ************************************************************************** !
!
! MatrixBlockAuxVarCopy: Copys an auxiliary object
! author: Glenn Hammond
! date: 03/04/2010
!
! ************************************************************************** !
subroutine MatrixBlockAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option  
  
  aux_var%dtotal = aux_var2%dtotal
  
end subroutine MatrixBlockAuxVarCopy

! ************************************************************************** !
!
! MatrixBlockAuxVarDestroy: Deallocates a matrix block auxiliary object
! author: Glenn Hammond
! date: 03/04/2010
!
! ************************************************************************** !
subroutine MatrixBlockAuxVarDestroy(aux_var)

  implicit none

  type(matrix_block_auxvar_type), pointer :: aux_var
  
  if (.not.associated(aux_var)) return
  
  if (associated(aux_var%dtotal))deallocate(aux_var%dtotal)
  nullify(aux_var%dtotal)
  
  deallocate(aux_var)
  nullify(aux_var)
  
end subroutine MatrixBlockAuxVarDestroy

! ************************************************************************** !
!
! MatrixBlockInfoCreate: Allocate and initialize matrix block info object
! author: Glenn Hammond
! date: 03/09/2010
!
! ************************************************************************** !
function MatrixBlockInfoCreate(dim1,dim2,dim3,option)

  use Option_module

  implicit none
  
  PetscInt :: dim1
  PetscInt :: dim2
  PetscInt :: dim3
  type(option_type) :: option
  
  type(matrix_block_info_type), pointer :: MatrixBlockInfoCreate
  
  type(matrix_block_info_type), pointer :: info

  allocate(info) 
  info%dim1 = dim1 
  info%dim2 = dim2 
  info%dim3 = dim3 

  MatrixBlockInfoCreate => info
  
end function MatrixBlockInfoCreate


! ************************************************************************** !
!
! MatrixBlockInfoDestroy: Deallocates a matrix block info object
! author: Glenn Hammond
! date: 03/08/2010
!
! ************************************************************************** !
subroutine MatrixBlockInfoDestroy(matrix_info)

  implicit none

  type(matrix_block_info_type), pointer :: matrix_info
  
  if (.not.associated(matrix_info)) return
  
  if (associated(matrix_info))deallocate(matrix_info)
  nullify(matrix_info)
  
end subroutine MatrixBlockInfoDestroy

end module Matrix_Block_Aux_module
