!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_out.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_out.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:35:51  lichtner
! Revised output formatting.
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

module ptran_out_module

  use ptran_global_module
  use trdynmem_module
  use fileio_module
  use ptran_speciation_module

private
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

public ptran_out, ptran_psi_out
contains
 
  
  

  subroutine ptran_out (kplt,da,da_1dof,da_kin)

  implicit none 


  DA  :: da, da_1dof, da_kin
  Vec :: c_nat, c_seq, por_nat, por_seq, pk_nat, pk_seq, rte_nat, rte_seq
  VecScatter :: scatter

  integer :: i,ibrk,icall,ierr,iflgam,ii,ix,j,jj,jy,k,kplt,kz,n,nn,fid,nr
  character(len=20) :: fname
  character(len=1) :: q
  real*8 :: xi, yj, zk
  real*8 :: sum1v, sum2v, crossarea, vol, phx !, fach2o, prod !, vel
! real*8 :: sum1, sum2
  real*8, pointer :: pk_p(:), rte_p(:)
  real*8, pointer :: fldvol(:)!, fldflx(:)
  real*8, pointer :: comp(:,:),cpsi(:,:),cmplx(:,:),alogjn(:), &
                     cc0(:),cx0(:),gam0(:),gamx0(:),dgamdi(:),ph(:)
  
  data icall/1/
  
  save icall

! -----------------------------------------
  
  if ((t < tplot(kplt) .and. ibrkcrv == 0) .or. iprint == -1) return
  
  if (ibrkcrv > 0 .and. icall == 1) then
    icall = 0
    if (myrank == 0) then
      write(fname,'(a9,a4)') 'ptran_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' kplot= ',kplt
      open(unit=ptran_IUNIT4,file=fname,action="write")
!     write(ptran_IUNIT4,'("#t            dt [",a1,"]",100i12)') tunit, &
!     (i,i=1,ibrkcrv), (i,i=1,ibrkcrv)
      if (jph == 0 .and. ihco3 == 0 .and. ncmplx == 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ((trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      else if (jph == 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ((trim(nam(j)),j=1,ncomp),(trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ('     pH     ',(trim(nam(j)),j=1,ncomp),(trim(nam(j)),j=1,ncomp), &
        i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 > 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ('     pH     ','    HCO3-   ',(trim(nam(j)),j=1,ncomp), &
        (trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      endif
    endif
  endif

  ierr = 0
  
  call DACreateNaturalVector(da,c_nat,ierr)
  call DAGlobalToNaturalBegin(da,c,INSERT_VALUES,c_nat,ierr)
  call DAGlobalToNaturalEnd(da,c,INSERT_VALUES,c_nat,ierr)

  call VecScatterCreateToAll(c_nat,scatter,c_seq,ierr)
  call VecScatterBegin(c_nat, c_seq, INSERT_VALUES, SCATTER_FORWARD, &
       scatter, ierr)
  call VecScatterEnd(c_nat, c_seq, INSERT_VALUES, SCATTER_FORWARD, &
       scatter, ierr)

! velocity fields
! call DACreateNaturalVector(da_3np, vl_nat, ierr)
! call DAGlobalToNaturalBegin(da_3np,vl,INSERT_VALUES,vl_nat,ierr)
! call DAGlobalToNaturalEnd(da_3np, vl, INSERT_VALUES,vl_nat,ierr)
! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr);  CHKERRQ(ierr)

  if (myrank == 0) call VecGetArrayF90(c_seq, c_p, ierr)

  if (ibrkcrv > 0) then
    if (myrank == 0) then
!     call VecGetArrayF90(vl_all, vl_p, ierr)
!     allocate(fldflx(ibrkcrv))
      allocate(ph(ibrkcrv))
      allocate(fldvol(ibrkcrv))
      allocate(comp(ncomp,ibrkcrv))
      allocate(cpsi(ncomp,ibrkcrv))
      allocate(cc0(ncomp))
      allocate(gam0(ncomp))
      allocate(dgamdi(ncomp))
      if (ncmplx > 0) then
        allocate(cmplx(ncmplx,ibrkcrv))
        allocate(alogjn(ncomp))
        allocate(cx0(ncmplx))
        allocate(gamx0(ncmplx))
      endif
      do ibrk = 1, ibrkcrv
!       fldflx(ibrk) = 0.d0
        fldvol(ibrk) = 0.d0
!       sum1 = 0.d0
!       sum2 = 0.d0
        sum1v = 0.d0
        sum2v = 0.d0
        do k = k1brk(ibrk),k2brk(ibrk)
          do j = j1brk(ibrk),j2brk(ibrk)
            do i = i1brk(ibrk),i2brk(ibrk)
              n = i+(j-1)*nx+(k-1)*nxy
              vol = dx0(i) * dy0(j) * dz0(k)

              if (ibrkface(ibrk) == 1) then
                nn = i-1+(j-1)*nx+(k-1)*nxy
                crossarea = dy0(j)*dz0(k)
              else if (ibrkface(ibrk) == 2) then
                nn = i+(j-2)*nx+(k-1)*nxy
                crossarea = dx0(i)*dz0(k)
              else if (ibrkface(ibrk) == 3) then
                nn = i+(j-1)*nx+(k-2)*nxy
                crossarea = dx0(i)*dy0(j)
              endif
              
!             vel = vl_p(ibrkface(ibrk)+3*(nn-1))

!             sum1 = sum1 + vel*crossarea*c_p(n)
!             sum2 = sum2 + vel*crossarea
              sum1v = sum1v + vol*c_p(1+(n-1)*ncomp) ! set to first primary species
              sum2v = sum2v + vol
!             print *,'output-brk: ',ibrk,i,j,k,n,ibrkface(ibrk),vol,c_p(n)
!             print *,'output-brk: ',ibrk,i,j,k,n,nn,ibrkface(ibrk), &
!             crossarea,vol,c_p(n),sum1v,sum2v !,vel*tconv,y_p(3*n), &
!             sum1,sum2
              do jj = 1, ncomp
                cc0(jj) = c_p(jj+(n-1)*ncomp)
                comp(jj,ibrk) = c_p(jj+(n-1)*ncomp)
                cpsi(jj,ibrk) = c_p(jj+(n-1)*ncomp)
                gam0(jj) = 1.d0
              enddo
              do ii = 1, ncmplx
                gamx0(ii) = 1.d0
              enddo

              !calculate activity coefficients and secondary species concentrations
              if (iact > 0) call trgameq(cc0,cx0,gam0,gamx0,dgamdi,iflgam)

!             print *,'ptran_out: ',adebye,bdebye
!             do jj = 1, ncomp
!               print *,'ptran_out: ',nam(jj),cc0(jj),gam0(jj),z(jj),a0(jj),bextend(jj)
!             enddo
!             do ii = 1, ncmplx
!               print *,'ptran_out: ',namcx(ii),cx0(ii),gamx0(ii),zx(ii),ax0(i), &
!               bextendx(ii)
!             enddo
              
              if(jph > 0) ph(ibrk) = -log10(c_p(jph+(n-1)*ncomp)*gam0(jph))

!             compute aqueous complex and total concentrations
              if (ncmplx > 0) then
                do ii = 1, ncmplx
                  cmplx(ii,ibrk) = cx0(ii)
                enddo
                do jj = 1, ncomp
                  do ii = 1, ncmplx
                    cpsi(jj,ibrk) = cpsi(jj,ibrk) + shom(jj,ii)*cmplx(ii,ibrk)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
!       if (sum2 .ne. 0.d0) fldflx(ibrk) = sum1/sum2
        if (sum2v .ne. 0.d0) fldvol(ibrk) = sum1v/sum2v
      enddo

      if (jph == 0 .and. jhco3 == 0 .and. ncmplx == 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        ((comp(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph == 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        ((comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        (ph(i),(comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 > 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        (ph(i),cmplx(ihco3,i),(comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp), &
        i=1,ibrkcrv)
      endif
      if (t < tplot(kplt)) call VecRestoreArrayF90(c_seq, c_p, ierr)
!     call VecRestoreArrayF90(vl_all, vl_p, ierr)
!     deallocate(fldflx)
      deallocate(fldvol)
    endif
  
    if (t < tplot(kplt)) then
      call VecDestroy(c_nat, ierr)
      call VecDestroy(c_seq, ierr)
!     call VecDestroy(vl_nat, ierr)
!     call VecDestroy(vl_all, ierr)
      return
    endif
  endif

  if (myrank==0) then
    q = ','

!---aqueous concentrations
    fid = 68
    if (kplt < 10) then
      write(fname,'(a4,i1,a4)') 'conc', kplt, '.dat'
    else
      write(fname,'(a4,i2,a4)') 'conc', kplt, '.dat'
    endif

    write(*,*) '--> write output file: ',fname

    open(unit=fid,file=fname,action="write")
    write(fid,10) t/tconv,tunit
    if (nx > 1 .and. ny == 1 .and. nz == 1) then
      if (jph == 0) then
        write(fid,21) (q,nam(j),j=1,ncomp)
      else
        write(fid,21) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,31) nx
      do ix=1,nx
        if (ix == 1) then
          xi = 0.5d0*dx0(ix)
        else
          xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
        endif
        n = ix
        if (jph == 0) then
          write(fid,40) xi,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        else
          phx = -log10(c_p(jph+(n-1)*ncomp))
          write(fid,40) xi,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        endif
      enddo    
    else if (nx == 1 .and. ny == 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,21) (q,nam(j),j=1,ncomp)
      else
        write(fid,21) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,31) nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        n = kz
        if (jph == 0) then
          write(fid,40) zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        else
          phx = -log10(c_p(jph+(n-1)*ncomp))
      
!       print *,'ptran_out: ',n,kplt,jph,zk,ph,log(10.d0),c_p(jph+(n-1)*ncomp)
      
          write(fid,40) zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        endif
      enddo    
    else if (nx > 1 .and. ny > 1 .and. nz == 1) then
      if (jph == 0) then
        write(fid,22) (q,nam(j),j=1,ncomp)
      else
        write(fid,22) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,32) nx,ny
      do jy=1,ny
        if (jy == 1) then
          yj = 0.5d0*dy0(jy)
        else
          yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix+(jy-1)*nx
          if (jph == 0) then
            write(fid,40) xi,yj,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          else
            phx = -log10(c_p(jph+(n-1)*ncomp))
            write(fid,40) xi,yj,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          endif
        enddo
      enddo
    else if (nx > 1 .and. ny == 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,23) (q,nam(j),j=1,ncomp)
      else
        write(fid,23) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,32) nx,nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix+(kz-1)*nx
          if (jph == 0) then
            write(fid,40) xi,zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          else
            phx = -log10(c_p(jph+(n-1)*ncomp))
            write(fid,40) xi,zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          endif
        enddo
      enddo
    else if (nx > 1 .and. ny > 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,24) (q,nam(j),j=1,ncomp)
      else
        write(fid,24) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,33) nx,ny,nz
      do kz = 1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do jy = 1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix = 1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx+(kz-1)*nxy
            if (jph == 0) then
              write(fid,40) xi,yj,zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
            else
              phx = -log10(c_p(jph+(n-1)*ncomp))
              write(fid,40) xi,yj,zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
            endif
          enddo
        enddo
      enddo
    endif
    close(fid)
    call VecRestoreArrayF90(c_seq,c_p,ierr)
  endif
  call VecDestroy(c_seq,ierr)
  call VecDestroy(c_nat,ierr)
    
  if (nkin > 0) then
  
!---mineral volume fractions & porosity
    call VecRestoreArrayF90(phik,phik_p,ierr) ! Is this needed-pcl?
    call VecRestoreArrayF90(por,por_p,ierr)

    call DACreateNaturalVector(da_kin,pk_nat,ierr)
    call DAGlobalToNaturalBegin(da_kin,phik,INSERT_VALUES,pk_nat,ierr)
    call DAGlobalToNaturalEnd(da_kin,phik,INSERT_VALUES,pk_nat,ierr)

    call VecScatterCreateToAll(pk_nat,scatter,pk_seq,ierr)
    call VecScatterBegin(pk_nat, pk_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(pk_nat, pk_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)

    call DACreateNaturalVector(da_1dof,por_nat,ierr)
    call DAGlobalToNaturalBegin(da_1dof,por,INSERT_VALUES,por_nat,ierr)
    call DAGlobalToNaturalEnd(da_1dof,por,INSERT_VALUES,por_nat,ierr)

    call VecScatterCreateToAll(por_nat,scatter,por_seq,ierr)
    call VecScatterBegin(por_nat, por_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(por_nat, por_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
  
    call VecGetArrayF90(por_seq,por_p,ierr)
    call VecGetArrayF90(pk_seq,pk_p,ierr)

    if (myrank == 0) then
      fid = 69
      if (kplt < 10) then
        write(fname,'(a4,i1,a4)') 'volf', kplt, '.dat'
      else
        write(fname,'(a4,i2,a4)') 'volf', kplt, '.dat'
      endif
      open(unit=fid,file=fname,action="write")
    
      if (myrank == 0) write(*,*) '--> write output file: ',fname
      
      write(fid,10) t/tconv,tunit

      if (nx > 1 .and. ny == 1 .and. nz == 1) then
        write(fid,21) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,31) nx
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix
          write(fid,40) xi,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
        enddo    
      else if (nx == 1 .and. ny == 1 .and. nz > 1) then
        write(fid,21) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,31) nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          n = kz
          write(fid,40) zk,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
        enddo    
      else if (nx > 1 .and. ny > 1 .and. nz == 1) then
        write(fid,22) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,32) nx,ny
        do jy=1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx
            write(fid,40) xi,yj,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
          enddo
        enddo
      else if (nx > 1 .and. ny == 1 .and. nz > 1) then
        write(fid,23) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,32) nx,nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(kz-1)*nx
            write(fid,40) xi,zk,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
          enddo
        enddo
      else if (nx > 1 .and. ny > 1 .and. nz > 1) then
        write(fid,24) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,33) nx,ny,nz
        do kz = 1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do jy = 1,ny
            if (jy == 1) then
              yj = 0.5d0*dy0(jy)
            else
              yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
            endif
            do ix = 1,nx
              if (ix == 1) then
                xi = 0.5d0*dx0(ix)
              else
                xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
              endif
              n = ix+(jy-1)*nx+(kz-1)*nxy
              write(fid,40) xi,yj,zk, &
              (pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
            enddo
          enddo
        enddo
      endif
      close(fid)
    endif
    call VecRestoreArrayF90(por_seq,por_p,ierr)
    call VecRestoreArrayF90(pk_seq,pk_p,ierr)
    call VecDestroy(pk_seq,ierr)
    call VecDestroy(pk_nat,ierr)
    call VecDestroy(por_seq,ierr)
    call VecDestroy(por_nat,ierr)

    call VecGetArrayF90(phik,phik_p,ierr)
    call VecGetArrayF90(por,por_p,ierr)

  
!---mineral reaction rates
    call VecRestoreArrayF90(rkin,rkin_p,ierr) ! Is this needed-pcl?

    call DACreateNaturalVector(da_kin,rte_nat,ierr)
    call DAGlobalToNaturalBegin(da_kin,rkin,INSERT_VALUES,rte_nat,ierr)
    call DAGlobalToNaturalEnd(da_kin,rkin,INSERT_VALUES,rte_nat,ierr)

    call VecScatterCreateToAll(rte_nat,scatter,rte_seq,ierr)
    call VecScatterBegin(rte_nat, rte_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(rte_nat, rte_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)

    call VecGetArrayF90(rte_seq,rte_p,ierr)

    if (myrank == 0) then
      fid = 69
      if (kplt < 10) then
        write(fname,'(a6,i1,a4)') 'rxnrte', kplt, '.dat'
      else
        write(fname,'(a6,i2,a4)') 'rxnrte', kplt, '.dat'
      endif
      open(unit=fid,file=fname,action="write")
    
      if (myrank == 0) write(*,*) '--> write output file: ',fname
      
      write(fid,10) t/tconv,tunit

      if (nx > 1 .and. ny == 1 .and. nz == 1) then
        write(fid,21) (q,namk(j),j=1,nkin)
        write(fid,31) nx
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix
          write(fid,40) xi,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
        enddo    
      else if (nx == 1 .and. ny == 1 .and. nz > 1) then
        write(fid,21) (q,namk(j),j=1,nkin)
        write(fid,31) nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          n = kz
          write(fid,40) zk,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
        enddo    
      else if (nx > 1 .and. ny > 1 .and. nz == 1) then
        write(fid,22) (q,namk(j),j=1,nkin)
        write(fid,32) nx,ny
        do jy=1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx
            write(fid,40) xi,yj,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
          enddo
        enddo
      else if (nx > 1 .and. ny == 1 .and. nz > 1) then
        write(fid,23) (q,namk(j),j=1,nkin)
        write(fid,32) nx,nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(kz-1)*nx
            write(fid,40) xi,zk,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
          enddo
        enddo
      else if (nx > 1 .and. ny > 1 .and. nz > 1) then
        write(fid,24) (q,namk(j),j=1,nkin)
        write(fid,33) nx,ny,nz
        do kz = 1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do jy = 1,ny
            if (jy == 1) then
              yj = 0.5d0*dy0(jy)
            else
              yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
            endif
            do ix = 1,nx
              if (ix == 1) then
                xi = 0.5d0*dx0(ix)
              else
                xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
              endif
              n = ix+(jy-1)*nx+(kz-1)*nxy
              write(fid,40) xi,yj,zk, &
              (rte_p(nr+(n-1)*nkin),nr=1,nkin)
            enddo
          enddo
        enddo
      endif
      close(fid)
    endif
    call VecRestoreArrayF90(rte_seq,rte_p,ierr)
    call VecDestroy(rte_seq,ierr)
    call VecDestroy(rte_nat,ierr)

    call VecGetArrayF90(rkin,rkin_p,ierr)
  endif

  10 format('TITLE= ','"Time= ',1pe12.4,' [',a1,']"')
  21 format('VARIABLES= " x/z "',100(a1,'"',a9'"'))
  22 format('VARIABLES= " x "," y "',100(a1,'"',a9'"'))
  23 format('VARIABLES= " x "," z "',100(a1,'"',a9'"'))
  24 format('VARIABLES= " x "," y "," z "',100(a1,'"',a9'"'))
  25 format('VARIABLES= " x "," y ","por"',100(a1,'"',a20'"'))
  30 format('ZONE T= " "',', I= ',i4,', J= ',i4,', K= ',i4)
  31 format('ZONE T= " "',', I= ',i4)
  32 format('ZONE T= " "',', I= ',i4,', J= ',i4)
  33 format('ZONE T= " "',', I= ',i4,', J= ',i4,', K= ',i4)
  40 format(1p100e12.4)

  end subroutine ptran_out
  
  
 !**********************************************************************************
 
   subroutine ptran_get_natural_psi(snpsi, ierr)
    implicit none
    
   real*8  :: snpsi(:)    
    integer :: ierr, ndex,na, status(MPI_STATUS_SIZE)
  integer :: ifound, jng, n
  real*8 pppsi(ncomp)
      
    ierr=-1
         ndex=1
    do na=0, nmax-1
        ifound=0
        if(myrank==0) then
         ifound=0
         do n=ndex,nlmax
         if(nL2A(n) == na ) then
          ifound=1
          pppsi= psi(:, nL2G(n))
          ndex=n
          exit   
        endif
         enddo
         if(ifound==0)then  
        call MPI_Recv(pppsi,ncomp, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, na, &
                PETSC_COMM_WORLD, status,ierr)           
         endif
        
        else
        do n=ndex,nlmax
          if(nL2A(n) == na) then
          jng = 1 + (nl2G(n)-1)*ncomp 
           call MPI_Send(psi(1:ncomp,nL2G(n)), ncomp ,MPI_DOUBLE_PRECISION, 0, na, &
                     PETSC_COMM_WORLD, ierr)   
            ndex=n
           exit
          endif
         enddo  
       endif            
    
      call MPI_Barrier(PETSC_COMM_WORLD, ierr)
      
      if(myrank ==0 )then
       snpsi(na*ncomp+1: (na+1)*ncomp) = pppsi(1:ncomp)
      endif
     

         
    enddo

   
   ierr=0
  
   end subroutine ptran_get_natural_psi    
  
 !**********************************************************************************
 
   subroutine ptran_get_natural_psig(snpsig, ierr)
    implicit none
    
   real*8  :: snpsig(:)    
    integer :: ierr, ndex,na, status(MPI_STATUS_SIZE)
  integer :: ifound, jng, n
  real*8 pppsig(ncomp)
      
    ierr=-1
         ndex=1
    do na=0, nmax-1
        ifound=0
        if(myrank==0) then
         ifound=0
         do n=ndex,nlmax
         if(nL2A(n) == na ) then
          ifound=1
          pppsig = psig(:, nL2G(n))
          ndex=n
          exit   
        endif
         enddo
         if(ifound==0)then  
        call MPI_Recv(pppsig, ncomp ,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, na+na, &
                PETSC_COMM_WORLD,status ,ierr) 
         endif
        
        else
        do n=ndex,nlmax
          if(nL2A(n) == na) then
          jng = 1 + (nl2G(n)-1)*ncomp 
             call MPI_Send(psig(1:ncomp,nL2G(n)), ncomp ,MPI_DOUBLE_PRECISION, 0, na+na, &
                     PETSC_COMM_WORLD, ierr)
            ndex=n
           exit
          endif
         enddo  
       endif            
    
      call MPI_Barrier(PETSC_COMM_WORLD, ierr)
      
      if(myrank ==0 )then
       snpsig(na*ncomp+1: (na+1)*ncomp) = pppsig(1:ncomp) 
      endif
     

         
    enddo

   
   ierr=0
  
   end subroutine ptran_get_natural_psig    

 !********************************************************************************** 
  subroutine ptran_psi_out (kplt,da,da_1dof,da_kin)

  implicit none 


  DA  :: da, da_1dof, da_kin
  Vec :: c_nat, c_seq, por_nat, por_seq, pk_nat, pk_seq, rte_nat, rte_seq
  VecScatter :: scatter

  integer :: i,ibrk,icall,ierr,iflgam,ii,ix,j,jj,jy,k,kplt,kz,n,nn,fid,nr
  character(len=20) :: fname
  character(len=1) :: q
  real*8 :: xi, yj, zk
  real*8 :: sum1v, sum2v, crossarea, vol, phx !, fach2o, prod !, vel
! real*8 :: sum1, sum2
  real*8, pointer :: pk_p(:), rte_p(:)
  real*8, pointer :: fldvol(:)!, fldflx(:)
  real*8, pointer :: comp(:,:),cpsi(:,:),cmplx(:,:),alogjn(:), &
                     cc0(:),cx0(:),gam0(:),gamx0(:),dgamdi(:),ph(:)
  real*8, allocatable :: snpsi(:), snpsig(:)  
  
  data icall/1/
  
  save icall

! -----------------------------------------
  
  if ((t < tplot(kplt) .and. ibrkcrv == 0) .or. iprint == -1) return
  
  if (ibrkcrv > 0 .and. icall == 1) then
    icall = 0
    if (myrank == 0) then
      write(fname,'(a9,a4)') 'ptran_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' kplot= ',kplt
      open(unit=ptran_IUNIT4,file=fname,action="write")
!     write(ptran_IUNIT4,'("#t            dt [",a1,"]",100i12)') tunit, &
!     (i,i=1,ibrkcrv), (i,i=1,ibrkcrv)
      if (jph == 0 .and. ihco3 == 0 .and. ncmplx == 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ((trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      else if (jph == 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ((trim(nam(j)),j=1,ncomp),(trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ('     pH     ',(trim(nam(j)),j=1,ncomp),(trim(nam(j)),j=1,ncomp), &
        i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 > 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'("#t            dt [",a1,"]",100a12)') tunit, &
        ('     pH     ','    HCO3-   ',(trim(nam(j)),j=1,ncomp), &
        (trim(nam(j)),j=1,ncomp),i=1,ibrkcrv)
      endif
    endif
  endif

  ierr = 0
  
  call DACreateNaturalVector(da,c_nat,ierr)
  call DAGlobalToNaturalBegin(da,c,INSERT_VALUES,c_nat,ierr)
  call DAGlobalToNaturalEnd(da,c,INSERT_VALUES,c_nat,ierr)

  call VecScatterCreateToAll(c_nat,scatter,c_seq,ierr)
  call VecScatterBegin(c_nat, c_seq, INSERT_VALUES, SCATTER_FORWARD, &
       scatter, ierr)
  call VecScatterEnd(c_nat, c_seq, INSERT_VALUES, SCATTER_FORWARD, &
       scatter, ierr)

! velocity fields
! call DACreateNaturalVector(da_3np, vl_nat, ierr)
! call DAGlobalToNaturalBegin(da_3np,vl,INSERT_VALUES,vl_nat,ierr)
! call DAGlobalToNaturalEnd(da_3np, vl, INSERT_VALUES,vl_nat,ierr)
! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr);  CHKERRQ(ierr)
 if (myrank == 0) then
   if (iphase >= 1) allocate(snpsi(nmax*ncomp))
   if (iphase  == 2) allocate(snpsig(nmax*ncomp))
 endif
  
  if (iphase >= 1) call ptran_get_natural_psi(snpsi, ierr)
  if (iphase == 2) call ptran_get_natural_psig(snpsig, ierr)
  
  if (myrank == 0) call VecGetArrayF90(c_seq, c_p, ierr)

  if (ibrkcrv > 0) then
    if (myrank == 0) then
!     call VecGetArrayF90(vl_all, vl_p, ierr)
!     allocate(fldflx(ibrkcrv))
      allocate(ph(ibrkcrv))
      allocate(fldvol(ibrkcrv))
      allocate(comp(ncomp,ibrkcrv))
      allocate(cpsi(ncomp,ibrkcrv))
      allocate(cc0(ncomp))
      allocate(gam0(ncomp))
      allocate(dgamdi(ncomp))
      if (ncmplx > 0) then
        allocate(cmplx(ncmplx,ibrkcrv))
        allocate(alogjn(ncomp))
        allocate(cx0(ncmplx))
        allocate(gamx0(ncmplx))
      endif
      do ibrk = 1, ibrkcrv
!       fldflx(ibrk) = 0.d0
        fldvol(ibrk) = 0.d0
!       sum1 = 0.d0
!       sum2 = 0.d0
        sum1v = 0.d0
        sum2v = 0.d0
        do k = k1brk(ibrk),k2brk(ibrk)
          do j = j1brk(ibrk),j2brk(ibrk)
            do i = i1brk(ibrk),i2brk(ibrk)
              n = i+(j-1)*nx+(k-1)*nxy
              vol = dx0(i) * dy0(j) * dz0(k)

              if (ibrkface(ibrk) == 1) then
                nn = i-1+(j-1)*nx+(k-1)*nxy
                crossarea = dy0(j)*dz0(k)
              else if (ibrkface(ibrk) == 2) then
                nn = i+(j-2)*nx+(k-1)*nxy
                crossarea = dx0(i)*dz0(k)
              else if (ibrkface(ibrk) == 3) then
                nn = i+(j-1)*nx+(k-2)*nxy
                crossarea = dx0(i)*dy0(j)
              endif
              
!             vel = vl_p(ibrkface(ibrk)+3*(nn-1))

!             sum1 = sum1 + vel*crossarea*c_p(n)
!             sum2 = sum2 + vel*crossarea
              sum1v = sum1v + vol*c_p(1+(n-1)*ncomp) ! set to first primary species
              sum2v = sum2v + vol
!             print *,'output-brk: ',ibrk,i,j,k,n,ibrkface(ibrk),vol,c_p(n)
!             print *,'output-brk: ',ibrk,i,j,k,n,nn,ibrkface(ibrk), &
!             crossarea,vol,c_p(n),sum1v,sum2v !,vel*tconv,y_p(3*n), &
!             sum1,sum2
              do jj = 1, ncomp
                cc0(jj) = c_p(jj+(n-1)*ncomp)
                comp(jj,ibrk) = c_p(jj+(n-1)*ncomp)
                cpsi(jj,ibrk) = c_p(jj+(n-1)*ncomp)
                gam0(jj) = 1.d0
              enddo
              do ii = 1, ncmplx
                gamx0(ii) = 1.d0
              enddo

              !calculate activity coefficients and secondary species concentrations
              if (iact > 0) call trgameq(cc0,cx0,gam0,gamx0,dgamdi,iflgam)

!             print *,'ptran_out: ',adebye,bdebye
!             do jj = 1, ncomp
!               print *,'ptran_out: ',nam(jj),cc0(jj),gam0(jj),z(jj),a0(jj),bextend(jj)
!             enddo
!             do ii = 1, ncmplx
!               print *,'ptran_out: ',namcx(ii),cx0(ii),gamx0(ii),zx(ii),ax0(i), &
!               bextendx(ii)
!             enddo
              
              if(jph > 0) ph(ibrk) = -log10(c_p(jph+(n-1)*ncomp)*gam0(jph))

!             compute aqueous complex and total concentrations
              if (ncmplx > 0) then
                do ii = 1, ncmplx
                  cmplx(ii,ibrk) = cx0(ii)
                enddo
                do jj = 1, ncomp
                  do ii = 1, ncmplx
                    cpsi(jj,ibrk) = cpsi(jj,ibrk) + shom(jj,ii)*cmplx(ii,ibrk)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
!       if (sum2 .ne. 0.d0) fldflx(ibrk) = sum1/sum2
        if (sum2v .ne. 0.d0) fldvol(ibrk) = sum1v/sum2v
      enddo

      if (jph == 0 .and. jhco3 == 0 .and. ncmplx == 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        ((comp(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph == 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        ((comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 == 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        (ph(i),(comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp),i=1,ibrkcrv)
      else if (jph > 0 .and. ihco3 > 0 .and. ncmplx > 0) then
        write(ptran_IUNIT4,'(1p100e12.4)') t/tconv,dt/tconv, &
!       (fldflx(i),i=1,ibrkcrv), &
!       (fldvol(i),i=1,ibrkcrv)
        (ph(i),cmplx(ihco3,i),(comp(j,i),j=1,ncomp),(cpsi(j,i),j=1,ncomp), &
        i=1,ibrkcrv)
      endif
      if (t < tplot(kplt)) call VecRestoreArrayF90(c_seq, c_p, ierr)
!     call VecRestoreArrayF90(vl_all, vl_p, ierr)
!     deallocate(fldflx)
      deallocate(fldvol)
    endif
  
    if (t < tplot(kplt)) then
      call VecDestroy(c_nat, ierr)
      call VecDestroy(c_seq, ierr)
!     call VecDestroy(vl_nat, ierr)
!     call VecDestroy(vl_all, ierr)
      return
    endif
  endif

  if (myrank==0) then
    q = ','

!---aqueous concentrations
    fid = 68
    if (kplt < 10) then
      write(fname,'(a4,i1,a4)') 'conc', kplt, '.dat'
    else
      write(fname,'(a4,i2,a4)') 'conc', kplt, '.dat'
    endif

    write(*,*) '--> write output file: ',fname

    open(unit=fid,file=fname,action="write")
    write(fid,10) t/tconv,tunit
    if (nx > 1 .and. ny == 1 .and. nz == 1) then
      if (jph == 0) then
        write(fid,21) (q,nam(j),j=1,ncomp)
      else
        write(fid,21) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,31) nx
      do ix=1,nx
        if (ix == 1) then
          xi = 0.5d0*dx0(ix)
        else
          xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
        endif
        n = ix
        if (jph == 0) then
          write(fid,40) xi,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        else
          phx = -log10(c_p(jph+(n-1)*ncomp))
          write(fid,40) xi,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        endif
      enddo    
    else if (nx == 1 .and. ny == 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,21) (q,nam(j),j=1,ncomp)
      else
        write(fid,21) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,31) nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        n = kz
        if (jph == 0) then
          write(fid,40) zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        else
          phx = -log10(c_p(jph+(n-1)*ncomp))
      
!       print *,'ptran_out: ',n,kplt,jph,zk,ph,log(10.d0),c_p(jph+(n-1)*ncomp)
      
          write(fid,40) zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
        endif
      enddo    
    else if (nx > 1 .and. ny > 1 .and. nz == 1) then
      if (jph == 0) then
        write(fid,22) (q,nam(j),j=1,ncomp)
      else
        write(fid,22) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,32) nx,ny
      do jy=1,ny
        if (jy == 1) then
          yj = 0.5d0*dy0(jy)
        else
          yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix+(jy-1)*nx
          if (jph == 0) then
            write(fid,40) xi,yj,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          else
            phx = -log10(c_p(jph+(n-1)*ncomp))
            write(fid,40) xi,yj,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          endif
        enddo
      enddo
    else if (nx > 1 .and. ny == 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,23) (q,nam(j),j=1,ncomp)
      else
        write(fid,23) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,32) nx,nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix+(kz-1)*nx
          if (jph == 0) then
            write(fid,40) xi,zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          else
            phx = -log10(c_p(jph+(n-1)*ncomp))
            write(fid,40) xi,zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
          endif
        enddo
      enddo
    else if (nx > 1 .and. ny > 1 .and. nz > 1) then
      if (jph == 0) then
        write(fid,24) (q,nam(j),j=1,ncomp)
      else
        write(fid,24) q,'pH',(q,nam(j),j=1,ncomp)
      endif
      write(fid,33) nx,ny,nz
      do kz = 1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do jy = 1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix = 1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx+(kz-1)*nxy
            if (jph == 0) then
              write(fid,40) xi,yj,zk,(c_p(j+(n-1)*ncomp),j=1,ncomp)
            else
              phx = -log10(c_p(jph+(n-1)*ncomp))
              write(fid,40) xi,yj,zk,phx,(c_p(j+(n-1)*ncomp),j=1,ncomp)
            endif
          enddo
        enddo
      enddo
    endif
    close(fid)
    call VecRestoreArrayF90(c_seq,c_p,ierr)
  endif
  call VecDestroy(c_seq,ierr)
  call VecDestroy(c_nat,ierr)
 
 
   ! write psi in aqeous phase
    if (myrank==0) then
    q = ','

  fid = 68
    if (kplt < 10) then
      write(fname,'(a4,i1,a4)') 'psi', kplt, '.dat'
    else
      write(fname,'(a4,i2,a4)') 'psi', kplt, '.dat'
    endif

    write(*,*) '--> write output file: ',fname

    open(unit=fid,file=fname,action="write")
    write(fid,10) t/tconv,tunit
    if (nx > 1 .and. ny == 1 .and. nz == 1) then
       write(fid,21) (q,nam(j),j=1,ncomp)
       write(fid,31) nx
      do ix=1,nx
        if (ix == 1) then
          xi = 0.5d0*dx0(ix)
        else
          xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
        endif
        n = ix
         write(fid,40) xi,(snpsi(j+(n-1)*ncomp),j=1,ncomp),&
     (snpsig(j+(n-1)*ncomp),j=1,ncomp)
         enddo    
     else if (nx == 1 .and. ny == 1 .and. nz > 1) then
        write(fid,21) (q,nam(j),j=1,ncomp)
         write(fid,31) nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        n = kz
         write(fid,40) zk,(snpsi(j+(n-1)*ncomp),j=1,ncomp)
      enddo    
    else if (nx > 1 .and. ny > 1 .and. nz == 1) then
       write(fid,22) (q,nam(j),j=1,ncomp)
       write(fid,32) nx,ny
      do jy=1,ny
        if (jy == 1) then
          yj = 0.5d0*dy0(jy)
        else
          yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix+(jy-1)*nx
          write(fid,40) xi,yj,(snpsi(j+(n-1)*ncomp),j=1,ncomp)
        enddo
      enddo
    else if (nx > 1 .and. ny == 1 .and. nz > 1) then
       write(fid,23) (q,nam(j),j=1,ncomp)
       write(fid,32) nx,nz
      do kz=1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
            n = ix+(kz-1)*nx
            write(fid,40) xi,zk,(snpsi(j+(n-1)*ncomp),j=1,ncomp)
         enddo
      enddo
    else if (nx > 1 .and. ny > 1 .and. nz > 1) then
         write(fid,24) (q,nam(j),j=1,ncomp)
         write(fid,33) nx,ny,nz
       do kz = 1,nz
        if (kz == 1) then
          zk = 0.5d0*dz0(kz)
        else
          zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
        endif
        do jy = 1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix = 1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx+(kz-1)*nxy
              write(fid,40) xi,yj,zk,(snpsi(j+(n-1)*ncomp),j=1,ncomp),&
        (snpsig(j+(n-1)*ncomp),j=1,ncomp)
            enddo
        enddo
      enddo
    endif
    close(fid)
  endif
 if (myrank == 0) then
   if (iphase >= 1) deallocate(snpsi)
   if (iphase == 2) deallocate(snpsig)
 endif
 

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (nkin > 0) then
  
!---mineral volume fractions & porosity
    call VecRestoreArrayF90(phik,phik_p,ierr) ! Is this needed-pcl?
    call VecRestoreArrayF90(por,por_p,ierr)

    call DACreateNaturalVector(da_kin,pk_nat,ierr)
    call DAGlobalToNaturalBegin(da_kin,phik,INSERT_VALUES,pk_nat,ierr)
    call DAGlobalToNaturalEnd(da_kin,phik,INSERT_VALUES,pk_nat,ierr)

    call VecScatterCreateToAll(pk_nat,scatter,pk_seq,ierr)
    call VecScatterBegin(pk_nat, pk_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(pk_nat, pk_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)

    call DACreateNaturalVector(da_1dof,por_nat,ierr)
    call DAGlobalToNaturalBegin(da_1dof,por,INSERT_VALUES,por_nat,ierr)
    call DAGlobalToNaturalEnd(da_1dof,por,INSERT_VALUES,por_nat,ierr)

    call VecScatterCreateToAll(por_nat,scatter,por_seq,ierr)
    call VecScatterBegin(por_nat, por_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(por_nat, por_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
  
    call VecGetArrayF90(por_seq,por_p,ierr)
    call VecGetArrayF90(pk_seq,pk_p,ierr)

    if (myrank == 0) then
      fid = 69
      if (kplt < 10) then
        write(fname,'(a4,i1,a4)') 'volf', kplt, '.dat'
      else
        write(fname,'(a4,i2,a4)') 'volf', kplt, '.dat'
      endif
      open(unit=fid,file=fname,action="write")
    
      if (myrank == 0) write(*,*) '--> write output file: ',fname
      
      write(fid,10) t/tconv,tunit

      if (nx > 1 .and. ny == 1 .and. nz == 1) then
        write(fid,21) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,31) nx
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix
          write(fid,40) xi,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
        enddo    
      else if (nx == 1 .and. ny == 1 .and. nz > 1) then
        write(fid,21) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,31) nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          n = kz
          write(fid,40) zk,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
        enddo    
      else if (nx > 1 .and. ny > 1 .and. nz == 1) then
        write(fid,22) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,32) nx,ny
        do jy=1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx
            write(fid,40) xi,yj,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
          enddo
        enddo
      else if (nx > 1 .and. ny == 1 .and. nz > 1) then
        write(fid,23) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,32) nx,nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(kz-1)*nx
            write(fid,40) xi,zk,(pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
          enddo
        enddo
      else if (nx > 1 .and. ny > 1 .and. nz > 1) then
        write(fid,24) (q,namk(j),j=1,nkin),q,'porosity'
        write(fid,33) nx,ny,nz
        do kz = 1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do jy = 1,ny
            if (jy == 1) then
              yj = 0.5d0*dy0(jy)
            else
              yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
            endif
            do ix = 1,nx
              if (ix == 1) then
                xi = 0.5d0*dx0(ix)
              else
                xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
              endif
              n = ix+(jy-1)*nx+(kz-1)*nxy
              write(fid,40) xi,yj,zk, &
              (pk_p(nr+(n-1)*nkin),nr=1,nkin),por_p(n)
            enddo
          enddo
        enddo
      endif
      close(fid)
    endif
    call VecRestoreArrayF90(por_seq,por_p,ierr)
    call VecRestoreArrayF90(pk_seq,pk_p,ierr)
    call VecDestroy(pk_seq,ierr)
    call VecDestroy(pk_nat,ierr)
    call VecDestroy(por_seq,ierr)
    call VecDestroy(por_nat,ierr)

    call VecGetArrayF90(phik,phik_p,ierr)
    call VecGetArrayF90(por,por_p,ierr)

  
!---mineral reaction rates
    call VecRestoreArrayF90(rkin,rkin_p,ierr) ! Is this needed-pcl?

    call DACreateNaturalVector(da_kin,rte_nat,ierr)
    call DAGlobalToNaturalBegin(da_kin,rkin,INSERT_VALUES,rte_nat,ierr)
    call DAGlobalToNaturalEnd(da_kin,rkin,INSERT_VALUES,rte_nat,ierr)

    call VecScatterCreateToAll(rte_nat,scatter,rte_seq,ierr)
    call VecScatterBegin(rte_nat, rte_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)
    call VecScatterEnd(rte_nat, rte_seq, INSERT_VALUES, SCATTER_FORWARD, &
         scatter, ierr)

    call VecGetArrayF90(rte_seq,rte_p,ierr)

    if (myrank == 0) then
      fid = 69
      if (kplt < 10) then
        write(fname,'(a6,i1,a4)') 'rxnrte', kplt, '.dat'
      else
        write(fname,'(a6,i2,a4)') 'rxnrte', kplt, '.dat'
      endif
      open(unit=fid,file=fname,action="write")
    
      if (myrank == 0) write(*,*) '--> write output file: ',fname
      
      write(fid,10) t/tconv,tunit

      if (nx > 1 .and. ny == 1 .and. nz == 1) then
        write(fid,21) (q,namk(j),j=1,nkin)
        write(fid,31) nx
        do ix=1,nx
          if (ix == 1) then
            xi = 0.5d0*dx0(ix)
          else
            xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
          endif
          n = ix
          write(fid,40) xi,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
        enddo    
      else if (nx == 1 .and. ny == 1 .and. nz > 1) then
        write(fid,21) (q,namk(j),j=1,nkin)
        write(fid,31) nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          n = kz
          write(fid,40) zk,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
        enddo    
      else if (nx > 1 .and. ny > 1 .and. nz == 1) then
        write(fid,22) (q,namk(j),j=1,nkin)
        write(fid,32) nx,ny
        do jy=1,ny
          if (jy == 1) then
            yj = 0.5d0*dy0(jy)
          else
            yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(jy-1)*nx
            write(fid,40) xi,yj,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
          enddo
        enddo
      else if (nx > 1 .and. ny == 1 .and. nz > 1) then
        write(fid,23) (q,namk(j),j=1,nkin)
        write(fid,32) nx,nz
        do kz=1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do ix=1,nx
            if (ix == 1) then
              xi = 0.5d0*dx0(ix)
            else
              xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
            endif
            n = ix+(kz-1)*nx
            write(fid,40) xi,zk,(rte_p(nr+(n-1)*nkin),nr=1,nkin)
          enddo
        enddo
      else if (nx > 1 .and. ny > 1 .and. nz > 1) then
        write(fid,24) (q,namk(j),j=1,nkin)
        write(fid,33) nx,ny,nz
        do kz = 1,nz
          if (kz == 1) then
            zk = 0.5d0*dz0(kz)
          else
            zk = zk+0.5d0*(dz0(kz)+dz0(kz-1))
          endif
          do jy = 1,ny
            if (jy == 1) then
              yj = 0.5d0*dy0(jy)
            else
              yj = yj+0.5d0*(dy0(jy)+dy0(jy-1))
            endif
            do ix = 1,nx
              if (ix == 1) then
                xi = 0.5d0*dx0(ix)
              else
                xi = xi+0.5d0*(dx0(ix)+dx0(ix-1))
              endif
              n = ix+(jy-1)*nx+(kz-1)*nxy
              write(fid,40) xi,yj,zk, &
              (rte_p(nr+(n-1)*nkin),nr=1,nkin)
            enddo
          enddo
        enddo
      endif
      close(fid)
    endif
    call VecRestoreArrayF90(rte_seq,rte_p,ierr)
    call VecDestroy(rte_seq,ierr)
    call VecDestroy(rte_nat,ierr)

    call VecGetArrayF90(rkin,rkin_p,ierr)
  endif

  10 format('TITLE= ','"Time= ',1pe12.4,' [',a1,']"')
  21 format('VARIABLES= " x/z "',100(a1,'"',a9'"'))
  22 format('VARIABLES= " x "," y "',100(a1,'"',a9'"'))
  23 format('VARIABLES= " x "," z "',100(a1,'"',a9'"'))
  24 format('VARIABLES= " x "," y "," z "',100(a1,'"',a9'"'))
  25 format('VARIABLES= " x "," y ","por"',100(a1,'"',a20'"'))
  30 format('ZONE T= " "',', I= ',i4,', J= ',i4,', K= ',i4)
  31 format('ZONE T= " "',', I= ',i4)
  32 format('ZONE T= " "',', I= ',i4,', J= ',i4)
  33 format('ZONE T= " "',', I= ',i4,', J= ',i4,', K= ',i4)
  40 format(1p100e12.4)

  end subroutine ptran_psi_out

  
end module ptran_out_module
