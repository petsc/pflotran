module rock_reac_module
  use pflow_gridtype_module
  use ptran_global_module
  use trdynmem_module
 
private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petsclog.h"
 
 public  Rock_reac
 
 contains
 subroutine Rock_reac(grid)
 
  use mixture_module  
   implicit none
    type(pflowGrid) :: grid 

 
  integer :: ierr
  integer*4 :: n
 


  PetscScalar, pointer :: yy_p(:),xx_p(:),iphase_p(:)
  
  PetscScalar, pointer ::  porosity_p(:), porosity0_p(:),volume_p(:) 
  PetscScalar, pointer :: perm_xx_p(:),  perm_yy_p(:),  perm_zz_p(:) 
  PetscScalar, pointer :: perm0_xx_p(:),  perm0_yy_p(:),  perm0_zz_p(:),&
                           perm_pow_p(:)          

!  integer, pointer ::iphase_p(:)
  
  
 
   !print *, 'updating grid'
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr)
  call VecGetArrayF90(grid%porosity0, porosity0_p, ierr)
  call VecGetArrayF90(grid%perm0_xx, perm0_xx_p, ierr)
  call VecGetArrayF90(grid%perm0_yy, perm0_yy_p, ierr)
  call VecGetArrayF90(grid%perm0_zz, perm0_zz_p, ierr)
  
  call VecGetArrayF90(grid%perm_pow, perm_pow_p, ierr)


  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
 call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  !call VecGetArrayF90(grid%icap, icap_p, ierr)
  !call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
 
  


 do n = 1, grid%nlmax
   !   if(iphase_p(n)==6 .or. iphase_p(n)==2)then
      perm_xx_p(n)=perm0_xx_p(n)*&
      (porosity_p(n)/porosity0_p(n))**perm_pow_p(n) 
        perm_yy_p(n)=perm0_yy_p(n)*&
      (porosity_p(n)/porosity0_p(n))**perm_pow_p(n) 
        perm_zz_p(n)=perm0_zz_p(n)*&
      (porosity_p(n)/porosity0_p(n))**perm_pow_p(n) 
    !  endif
 enddo
      
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr)
  call VecRestoreArrayF90(grid%porosity0, porosity_p, ierr)
  call VecRestoreArrayF90(grid%perm0_xx, perm0_xx_p, ierr)
  call VecRestoreArrayF90(grid%perm0_yy, perm0_yy_p, ierr)
  call VecRestoreArrayF90(grid%perm0_zz, perm0_zz_p, ierr)
 
  call VecRestoreArrayF90(grid%perm_pow, perm_pow_p, ierr)


  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)

 ! call VecRestoreArrayF90(grid%icap, icap_p, ierr)
 ! call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)


 
end subroutine Rock_reac
end module rock_reac_module
