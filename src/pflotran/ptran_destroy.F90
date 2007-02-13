!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_destroy.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_destroy.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:29:50  lichtner
! Removed PetscFinalize to main program.
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module ptran_destroy_module

contains

  subroutine ptran_destroy (da,da_mat,da_1dof,da_kin)
   
  use ptran_global_module
  use trdynmem_module
   
  implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
!#include "include/finclude/petscviewer.h"

  DA    :: da ,da_mat, da_1dof, da_kin
  
  integer :: ierr
  
  ierr = 0

  call DADestroy(da,ierr)
  call DADestroy(da_mat,ierr)
  call DADestroy(da_1dof,ierr)

  call MatDestroy(A,ierr)
  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  call VecDestroy(c,ierr)
  call VecDestroy(cc,ierr)
  call VecDestroy(dc,ierr)
  call VecDestroy(ccloc,ierr)
  call VecDestroy(porosity,ierr)
  call VecDestroy(porloc,ierr)
  call VecDestroy(tortuosity,ierr)
  call VecDestroy(tort_loc,ierr)
  call VecDestroy(temp,ierr)
  call VecDestroy(temploc,ierr)
  call VecDestroy(press,ierr)
  call VecDestroy(pressloc,ierr)
  
  call VecRestoreArrayF90(ghost_loc,ghost_loc_p,ierr)
  call VecDestroy(ghost_loc,ierr)

  if (nkin > 0) then
    call VecRestoreArrayF90(phik,phik_p,ierr)
    call VecRestoreArrayF90(surf,surf_p,ierr)
    call VecRestoreArrayF90(phik0,phik0_p,ierr)
    call VecRestoreArrayF90(surf0,surf0_p,ierr)
    call VecRestoreArrayF90(por,por_p,ierr)
    call VecRestoreArrayF90(rkin,rkin_p,ierr)
    call VecRestoreArrayF90(rrkin,rrkin_p,ierr)
    call DADestroy(da_kin,ierr)
    call VecDestroy(phik,ierr)
    call VecDestroy(surf,ierr)
    call VecDestroy(phik0,ierr)
    call VecDestroy(surf0,ierr)
    call VecDestroy(rkin,ierr)
    call VecDestroy(rrkin,ierr)
  endif

! call VecDestroy(siteden,ierr)
! call VecDestroy(csorpf,ierr)
! call VecDestroy(csorp,ierr)
! call VecDestroy(ccsorp,ierr)
    
  end subroutine ptran_destroy

end module ptran_destroy_module
