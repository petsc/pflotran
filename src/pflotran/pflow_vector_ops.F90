! pflow_vector_ops.F90
! This file contains functions/subroutines for manipulating PETSc vectors.

! pflow_pack_xx()
! This subroutine fills the vector xx in the following manner:
! The first v1_dof elements of v1, followed by the first v2_dof elements
! of v2, then the next v1_dof elements of v1, and so on.
! Eventually this subroutine should be extended to handle any number of 
! vectors from which xx is to be filled.
!
! NOTE: I believe xx, v1, and v2 must already have been assembled. 

module pflow_vector_ops_module

  implicit none

  private  
#include "definitions.h"

  public :: pflow_pack_xx2, pflow_pack_xx3, pflow_pack_xx4
contains

  subroutine pflow_pack_xx2(xx, v1, v1_dof, v2, v2_dof, ierr)
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec, intent(inout) :: xx
  Vec, intent(in), optional :: v1, v2
  PetscInt, intent(in), optional :: v1_dof, v2_dof
  PetscErrorCode, intent(out), optional :: ierr

  PetscInt :: i
  PetscReal, pointer :: xx_p(:), v1_p(:), v2_p(:)
  PetscInt :: xx_i, v1_i, v2_i  ! Used to index entries of xx, v1, v2.
  PetscInt :: xx_size, v1_size, v2_size

  call VecGetLocalSize(xx, xx_size, ierr)
  call VecGetLocalSize(v1, v1_size, ierr)
  call VecGetLocalSize(v2, v2_size, ierr)
  if(xx_size .ne. v1_size + v2_size) then
    print *, "Error in pflow_pack_xx: Vector sizes are inconsistent."
    ierr = 1
    return
  endif

  ! Better put some error checking here eventually (make sure vectors are
  ! the right sizes).
  
  call VecGetArrayF90(xx, xx_p, ierr)
  call VecGetArrayF90(v1, v1_p, ierr)
  call VecGetArrayF90(v2, v2_p, ierr)

  xx_i = 0; v1_i = 0; v2_i = 0
  do ! while(xx_i < xx_size)
    do i=1, v1_dof
      v1_i = v1_i + 1
      xx_i = xx_i + 1
      xx_p(xx_i) = v1_p(v1_i)
    enddo

    do i=1, v2_dof
      v2_i = v2_i + 1
      xx_i = xx_i + 1
      xx_p(xx_i) = v2_p(v2_i)
    enddo
    if (xx_i == xx_size) exit
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr)
  call VecRestoreArrayF90(v1, v1_p, ierr)
  call VecRestoreArrayF90(v2, v2_p, ierr)

  end subroutine pflow_pack_xx2

  subroutine pflow_pack_xx3(xx, v1, v1_dof, v2, v2_dof, v3, v3_dof, ierr)
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  Vec, intent(inout) :: xx
  Vec, intent(in) :: v1, v2, v3
  PetscInt, intent(in) :: v1_dof, v2_dof, v3_dof
  PetscErrorCode, intent(out) :: ierr

! Vec, intent(in), optional :: v1, v2, v3
! PetscInt, intent(in), optional :: v1_dof, v2_dof, v3_dof
! PetscInt, intent(out), optional :: ierr

  PetscInt :: i
  PetscReal, pointer :: xx_p(:), v1_p(:), v2_p(:), v3_p(:)
  PetscInt :: xx_i, v1_i, v2_i, v3_i  ! Used to index entries of xx, v1, v2, v3
  PetscInt :: xx_size, v1_size, v2_size, v3_size

! call VecView(v3,PETSC_VIEWER_STDOUT_WORLD,ierr)

  call VecGetLocalSize(xx, xx_size, ierr)
  call VecGetLocalSize(v1, v1_size, ierr)
  call VecGetLocalSize(v2, v2_size, ierr)
  call VecGetLocalSize(v3, v3_size, ierr)
  
! print *,'pflow_vector_ops: ',v1_dof,v2_dof,v3_dof,v1_size,v2_size,v3_size
  
  if(xx_size .ne. v1_size + v2_size + v3_size) then
    print *, "Error in pflow_pack_xx: Vector sizes are inconsistent.xx3"
    ierr = 1
    return
  endif

  ! Better put some error checking here eventually (make sure vectors are
  ! the right sizes).
  
  call VecGetArrayF90(xx, xx_p, ierr)
  call VecGetArrayF90(v1, v1_p, ierr)
  call VecGetArrayF90(v2, v2_p, ierr)
  call VecGetArrayF90(v3, v3_p, ierr)

  xx_i = 0; v1_i = 0; v2_i = 0; v3_i = 0
  do ! while(xx_i < xx_size)
    do i=1, v1_dof
      v1_i = v1_i + 1
      xx_i = xx_i + 1
      xx_p(xx_i) = v1_p(v1_i)
    enddo

    do i=1, v2_dof
      v2_i = v2_i + 1
      xx_i = xx_i + 1
      xx_p(xx_i) = v2_p(v2_i)
    enddo

    do i=1, v3_dof
      v3_i = v3_i + 1
      xx_i = xx_i + 1
      xx_p(xx_i) = v3_p(v3_i)
    enddo
    if (xx_i == xx_size) exit
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr)
  call VecRestoreArrayF90(v1, v1_p, ierr)
  call VecRestoreArrayF90(v2, v2_p, ierr)
  call VecRestoreArrayF90(v3, v3_p, ierr)

  end subroutine pflow_pack_xx3




 subroutine pflow_pack_xx4(xx, v1, v1_dof, v2, v2_dof, v3, v3_dof,&
                           v4, v4_dof, ierr)
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  Vec, intent(inout) :: xx
  Vec, intent(in) :: v1, v2, v3, v4
  PetscInt, intent(in) :: v1_dof, v2_dof, v3_dof,v4_dof
  PetscErrorCode, intent(out) :: ierr

! Vec, intent(in), optional :: v1, v2, v3
! PetscInt, intent(in), optional :: v1_dof, v2_dof, v3_dof
! PetscInt, intent(out), optional :: ierr

  PetscReal, pointer :: xx_p(:), v1_p(:), v2_p(:), v3_p(:),v4_p(:)
  PetscInt :: xx_i, v1_i, v2_i, v3_i, v4_i
  ! Used to index entries of xx, v1, v2, v3
  PetscInt :: xx_size, v1_size, v2_size, v3_size, v4_size

! call VecView(v3,PETSC_VIEWER_STDOUT_WORLD,ierr)

  call VecGetLocalSize(xx, xx_size, ierr)
  call VecGetLocalSize(v1, v1_size, ierr)
  call VecGetLocalSize(v2, v2_size, ierr)
  call VecGetLocalSize(v3, v3_size, ierr)
  call VecGetLocalSize(v4, v4_size, ierr)

! print *,'pflow_vector_ops: ',v1_dof,v2_dof,v3_dof,v1_size,v2_size,v3_size
  
!  if(xx_size .ne. v1_size + v2_size + v3_size + v4_size) then
!    print *, "Error in pflow_pack_xx: Vector sizes are inconsistent."
!    ierr = 1
!    return
!  endif

  ! Better put some error checking here eventually (make sure vectors are
  ! the right sizes).
  
  call VecGetArrayF90(xx, xx_p, ierr)
  call VecGetArrayF90(v1, v1_p, ierr)
  call VecGetArrayF90(v2, v2_p, ierr)
  call VecGetArrayF90(v3, v3_p, ierr)
  call VecGetArrayF90(v4, v4_p, ierr)

  xx_i = 0; v1_i = 0; v2_i = 0; v3_i = 0; v4_i = 0
  do ! while(xx_i < xx_size)
      v1_i = v1_i + v1_dof
      xx_i = xx_i + 1
   
      xx_p(xx_i) = v1_p(v1_i)

      v2_i = v2_i +  v2_dof
      xx_i = xx_i + 1
      xx_p(xx_i) = v2_p(v2_i)

      v3_i = v3_i + v3_dof
      xx_i = xx_i + 1
      xx_p(xx_i) = v3_p(v3_i)

      v4_i = v4_i + v4_dof
      xx_i = xx_i + 1
      xx_p(xx_i) = v4_p(v4_i)

!   print *,'pflow_vector_ops/pack4: ', v1_i/v1_dof, v1_p(v1_i), &
!   v2_p(v2_i),v3_p(v3_i), v4_p(v4_i)

    if (xx_i == xx_size) exit
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr)
  call VecRestoreArrayF90(v1, v1_p, ierr)
  call VecRestoreArrayF90(v2, v2_p, ierr)
  call VecRestoreArrayF90(v3, v3_p, ierr)
  call VecRestoreArrayF90(v4, v4_p, ierr)

end subroutine pflow_pack_xx4


  subroutine pflow_update_fldvar (xx, press, temp, sat, xmol, den,  &
           porosity, enthalpy, accum, dencpr, ithrm, iphas, scale, eqkair, &
           nlmax, ndof, nph, jgas, jh2o)

  use water_eos_module

  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec, intent(in) :: xx
  Vec, intent(inout) :: press, temp, sat, xmol, den, porosity, enthalpy, &
                        accum, ithrm, iphas
  PetscInt, intent(in) :: nlmax, ndof, nph, jgas, jh2o
  PetscReal, pointer :: ithrm_p(:), iphas_p(:), xx_p(:), &
                     sat_p(:), xmol_p(:), den_p(:), por_p(:), &
                     press_p(:), temp_p(:), enthalpy_p(:), accum_p(:)
  PetscInt :: i, j1, j2, j3, m
  PetscErrorCode :: ierr
  PetscReal :: dl,dg,sl,sg,pa,pg,pl,pv,ps,hl,hg,por,porm1,scale,tc, &
  eqkair,xl1,xl2,xg1,xg2,dencp,dencpr(*)

  call VecGetArrayF90(xx,    xx_p,    ierr)
  call VecGetArrayF90(press, press_p, ierr)
  call VecGetArrayF90(temp,  temp_p,  ierr)
  call VecGetArrayF90(sat,   sat_p,   ierr)
  call VecGetArrayF90(xmol,  xmol_p,  ierr)
  call VecGetArrayF90(den,   den_p,   ierr)
  call VecGetArrayF90(porosity, por_p, ierr)
  call VecGetArrayF90(enthalpy, enthalpy_p, ierr)
  call VecGetArrayF90(accum, accum_p, ierr)
  call VecGetArrayF90(ithrm, ithrm_p, ierr)
  call VecGetArrayF90(iphas, iphas_p, ierr)
  
  do m = 1, nlmax
    j1 = 1+(m-1)*ndof
    j2 = 2+(m-1)*ndof
    j3 = 3+(m-1)*ndof
    
    por = por_p(m)
    porm1 = 1.d0 - por
    
    i = ithrm_p(m)
    dencp = dencpr(i)

!   print *,'fldvar: ',m,iphas_p(m),por,sat_p(m)
    
    if (iphas_p(m) == 0) then ! pure gas phase

      ! unpack primary variables
      pg  = xx_p(j1)
      tc  = xx_p(j2)
      xg2 = xx_p(j3)

      ! store field variables
      press_p(jgas+(m-1)*nph) = pg
      temp_p(m)               = tc
      xmol_p(jgas+(m-1)*nph)  = xg2
      
      ! compute accumulation term at current time step
      dg = den_p(jgas+(m-1)*nph)
      hg = enthalpy_p(jgas+(m-1)*nph)
      xg1 = 1.d0 - xg2
      accum_p(j1) = por * dg * xg1
      accum_p(j2) = por * (dg * hg - scale * pg) + porm1 * dencp * tc
      accum_p(j3) = por * dg * xg2

      ! check for phase change
      pa = xg2 * pg
      pv = pg - pa
      call psat(tc,ps,ierr)
      if (pv >= ps) then
        iphas_p(m) = 2 ! gas -> 2ph
      endif

    else if (iphas_p(m) == 1) then ! pure liquid phase

      ! unpack primary variables
      pl  = xx_p(1+(m-1)*ndof)
      tc  = xx_p(2+(m-1)*ndof)
      xl2 = xx_p(3+(m-1)*ndof)
      
      ! store field variables
      press_p(jh2o+(m-1)*nph) = pl
      temp_p(m)               = tc
      xmol_p(jh2o+(m-1)*nph)  = xl2

      ! compute accumulation term at current time step
      dl = den_p(jh2o+(m-1)*nph)
      hl = enthalpy_p(jh2o+(m-1)*nph)
      xl1 = 1.d0 - xl2
      accum_p(j1) = por * dl * xl1
      accum_p(j2) = por * (dl * hl - scale * pl) + porm1 * dencp * tc
      accum_p(j3) = por * dl * xl2

      ! check for phase change
      pa = eqkair * xl2
      pv = pl - pa
      call psat(tc,ps,ierr)
      
!     if (pv <= ps) then
!       iphas_p(m) = 2 ! liq -> 2ph
!     endif

!     print *,'fldvar: ',m,iphas_p(m),tc,pl !,pa,xl2,eqkair,pv,ps

    else if (iphas_p(m) == 2) then ! 2 phase

      ! unpack primary variables
      pg  = xx_p(1+(m-1)*ndof)
      sl  = xx_p(2+(m-1)*ndof)
      xg2 = xx_p(3+(m-1)*ndof)

      ! store field variables
      press_p(jgas+(m-1)*nph) = pg
      sat_p(jh2o+(m-1)*nph)   = sl
      xmol_p(jgas+(m-1)*nph)  = xg2
      
      ! compute accumulation term at current time step
      xg1 = 1.d0 - xg2
      xl2 = pa / eqkair
      xl1 = 1.d0 - xl2
      sg = 1.d0 - sl
      pa = xg2 * pg
      pv = pg - pa
!     tc = tsat(pv)
      dl = den_p(jh2o+(m-1)*nph)
      dg = den_p(jgas+(m-1)*nph)
      hl = enthalpy_p(jh2o+(m-1)*nph)
      hg = enthalpy_p(jgas+(m-1)*nph)
!     pc = pcap(sl)
!     pl = pg - pc
      accum_p(j1) = por * (sl * dl * xl1 + sg * dg * xg1)
      accum_p(j2) = por * (sl * (dl * hl - scale * pl) &
        + sg * (dg * hg - scale * pg)) + porm1 * dencp * tc
      accum_p(j3) = por * (sl * dl * xl2 + sg * dg * xg2)

      ! check for phase change
      if (sat_p(jgas+(m-1)*nph) >= 1.d0) then
        iphas_p(m) = 0 ! 2ph -> gas
      else if (sat_p(jgas+(m-1)*nph) <= 0.d0) then
        iphas_p(m) = 1 ! 2ph -> liq
      endif
    endif
  enddo
  call VecRestoreArrayF90(xx,    xx_p,    ierr)
  call VecRestoreArrayF90(press, press_p, ierr)
  call VecRestoreArrayF90(temp,  temp_p,  ierr)
  call VecRestoreArrayF90(sat,   sat_p,   ierr)
  call VecRestoreArrayF90(xmol,  xmol_p,  ierr)
  call VecRestoreArrayF90(den,   den_p,   ierr)
  call VecRestoreArrayF90(porosity, por_p, ierr)
  call VecRestoreArrayF90(enthalpy, enthalpy_p, ierr)
  call VecRestoreArrayF90(accum, accum_p, ierr)
  call VecRestoreArrayF90(ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(iphas, iphas_p, ierr)

  end subroutine pflow_update_fldvar

end module pflow_vector_ops_module
