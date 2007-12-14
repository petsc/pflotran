!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_conn.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_conn.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:29:00  lichtner
! Minor change.
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

module ptran_conn_module

  public

  interface ptran_read_vel
    module procedure ptran_read_vel
  end interface ptran_read_vel

  interface ptran_conn
    module procedure ptran_conn
  end interface ptran_conn

contains

  subroutine ptran_conn (da_1dof)

  use ptran_global_module
  use trdynmem_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA :: da_1dof
  
  integer :: i,ierr,ir,ird,j,k,m,mg1,mg2,nc,ng
  integer :: ii1,ii2,jj1,jj2,kk1,kk2
  
  real*8, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  
  nconn = (ngx-1)*nly*nlz + nlx*(ngy-1)*nlz + nlx*nly*(ngz-1)

  allocate(vl(nconn))
  if (iphase == 2) allocate(vg(nconn))

  ! Extract local, ghosted portions of dx, dy, dz vectors.
  call DAGlobalToLocalBegin(da_1dof, dx, INSERT_VALUES, dx_loc, ierr)
  call DAGlobalToLocalEnd(da_1dof, dx, INSERT_VALUES, dx_loc, ierr)
  call DAGlobalToLocalBegin(da_1dof, dy, INSERT_VALUES, dy_loc, ierr)
  call DAGlobalToLocalEnd(da_1dof, dy, INSERT_VALUES, dy_loc,ierr)
  call DAGlobalToLocalBegin(da_1dof, dz, INSERT_VALUES, dz_loc, ierr)
  call DAGlobalToLocalEnd(da_1dof, dz, INSERT_VALUES, dz_loc,ierr)
  call VecGetArrayF90(dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(dz_loc, dz_loc_p, ierr)

!   connections are local/ghosted
!   interior nodes

!   if (using_pflowGrid .ne. PETSC_TRUE) then
    
      allocate(nd1(nconn))
      allocate(nd2(nconn))
      allocate(dist1(nconn))
      allocate(dist2(nconn))
      allocate(area(nconn))

      nc = 0

!-----x-connections
      if(ngx > 1) then
        do k = kstart,kend
          do j = jstart,jend
            do i = 1,ngx-1
              mg1 = i+j*ngx+k*ngxy
              mg2 = mg1+1
              nc = nc+1
              nd1(nc) = mg1
              nd2(nc) = mg2
              dist1(nc) = 0.5d0*dx_loc_p(mg1)
              dist2(nc) = 0.5d0*dx_loc_p(mg2)
              if (igeom == 1) then ! cartesian
                area(nc) = dy_loc_p(mg1)*dz_loc_p(mg1)
              else if (igeom == 2) then ! cylindrical
                area(nc) = 2.d0*Pi*rd(i+nxs)*dz_loc_p(mg1)
                !print *,'ptran_conn-x: ',myrank,nc,i,nxs,mg1,mg2, &
                !rd(i+nxs),dz_loc_p(mg1),area(nc)
              else if (igeom == 3) then ! spherical
                area(nc) = 4.d0*Pi*rd(i+nxs)**2
              endif
              vl(nc) = vlx0
              if (iphase == 2) vg(nc) = vgx0
            enddo
          enddo
        enddo
      endif
      
      nconnx = nc

!-----y-connections
      if(ngy > 1) then
        do k = kstart,kend
          do i = istart,iend
            do j = 1,ngy-1
              mg1 = i+1+(j-1)*ngx+k*ngxy
              mg2 = mg1+ngx
              nc = nc+1
              nd1(nc) = mg1
              nd2(nc) = mg2
              dist1(nc) = 0.5d0*dy_loc_p(mg1)
              dist2(nc) = 0.5d0*dy_loc_p(mg2)
              if (igeom == 1) then ! cartesian
                area(nc) = dx_loc_p(mg1)*dz_loc_p(mg1)
              endif
              vl(nc) = vly0
              if (iphase == 2) vg(nc) = vgy0
            enddo
          enddo
        enddo
      endif

      nconny = nc
      
!-----z-connections
      if(ngz > 1) then
        do j = jstart,jend
          do i = istart,iend
            do k = 1,ngz-1
              mg1 = i+1+j*ngx+(k-1)*ngxy
              mg2 = mg1+ngxy          
              nc = nc+1
              nd1(nc) = mg1
              nd2(nc) = mg2
              dist1(nc) = 0.5d0*dz_loc_p(mg1)
              dist2(nc) = 0.5d0*dz_loc_p(mg2)
              if (igeom == 1) then ! cartesian
                area(nc) = dx_loc_p(mg1)*dy_loc_p(mg1)
              else if (igeom == 2) then ! cylindrical
                area(nc) = Pi*(rd(i+nxs+1)+rd(i+nxs))* &
                (rd(i+nxs+1)-rd(i+nxs))
                !print *,'ptran_conn-z: ',myrank,nc,i,nxs,mg1,mg2, &
                !rd(i+nxs+1),rd(i+nxs),area(nc)
              endif
              vl(nc) = vlz0
              if (iphase == 2) vg(nc) = vgz0
            enddo
          enddo
        enddo
      endif

      if (nconn .ne. nc) then
        write(*,*) 'error in computing connection geometry: ', &
        'rank = ',myrank,' nconn = ',nconn,' nc =',nc
        stop
      endif

!   endif
    
!   boundary nodes: local\non-ghosted
    
    if (nx > 1 .and. ny == 1 .and. nz == 1) then
      if (nxs == ngxs) nconnbc = nconnbc + 1
      if (nxe == ngxe) nconnbc = nconnbc + 1
    else if (nx == 1 .and. ny == 1 .and. nz > 1) then
      if (nzs == ngzs) nconnbc = nconnbc + 1
      if (nze == ngze) nconnbc = nconnbc + 1
    else
      if (nx > 1) then
        if (nxs == ngxs) nconnbc = nconnbc + nlyz
        if (nxe == ngxe) nconnbc = nconnbc + nlyz
      endif
      if (ny > 1) then
        if (nys == ngys) nconnbc = nconnbc + nlxz
        if (nye == ngye) nconnbc = nconnbc + nlxz
      endif
      if (nz > 1) then
        if (nzs == ngzs) nconnbc = nconnbc + nlxy
        if (nze == ngze) nconnbc = nconnbc + nlxy
      endif
    endif
    if (nconnbc > 0) then
      allocate(mblkbc(nconnbc), ibconn(nconnbc))
!     if (using_pflowGrid .ne. PETSC_TRUE) then
        allocate(distbc(nconnbc))
        allocate(areabc(nconnbc))
!    allocate(xphibc(nconnbc))
!     endif
      allocate(vlbc(nconnbc))
      if (iphase == 2) allocate(vgbc(nconnbc))
    endif
    
    write(*,'(" --> ptranconn: rank = ",i4, &
   &", boundary connections =", i6)') myrank,nconnbc

    nc = 0
    
    if (nxs == ngxs .or. nxe == ngxe &
   .or. nys == ngys .or. nye == ngye &
   .or. nzs == ngzs .or. nze == ngze) then
   
      do ibc = 1, nblkbc
        do ir = iregbc1(ibc),iregbc2(ibc)

          kk1 = k1bc(ir) - nzs
          kk2 = k2bc(ir) - nzs
          jj1 = j1bc(ir) - nys
          jj2 = j2bc(ir) - nys
          ii1 = i1bc(ir) - nxs
          ii2 = i2bc(ir) - nxs

          kk1 = max(1,kk1)
          kk2 = min(nlz,kk2)
          jj1 = max(1,jj1)
          jj2 = min(nly,jj2)
          ii1 = max(1,ii1)
          ii2 = min(nlx,ii2)

          if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle

          do k = kk1,kk2
            do j = jj1,jj2
              do i = ii1,ii2
                nc = nc + 1
                m = i+(j-1)*nlx+(k-1)*nlxy
                mblkbc(nc) = m
                ibconn(nc) = ibc
                ng = nL2G(m)
                
                select case (igeom)
                
                case(1) ! cartesian
                
                if (iface(ibc) == 1) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                  vlbc(nc)   = vlx0
                  if (iphase == 2) vgbc(nc) = vgx0
                else if (iface(ibc) == 2) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                  vlbc(nc)   = -vlx0
                  if (iphase == 2) vgbc(nc) = -vgx0
                else if (iface(ibc) == 3) then
                  distbc(nc) = 0.5d0*dz_loc_p(ng)
                  areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                  vlbc(nc)   = vlz0
                  if (iphase == 2) vgbc(nc) = vgz0
                else if (iface(ibc) == 4) then
                  distbc(nc) = 0.5d0*dz_loc_p(ng)
                  areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                  vlbc(nc)   = -vlz0
                  if (iphase == 2) vgbc(nc) = -vgz0
                else if (iface(ibc) == 5) then
                  distbc(nc) = 0.5d0*dy_loc_p(ng)
                  areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                  vlbc(nc)   = vly0
                  if (iphase == 2) vgbc(nc) = vgy0
                else if (iface(ibc) == 6) then
                  distbc(nc) = 0.5d0*dy_loc_p(ng)
                  areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                  vlbc(nc)   = -vly0
                  if (iphase == 2) vgbc(nc) = -vgz0
                endif
                
                case (2) ! cylindrical
                
!               ird = mod(mod((m),nlxy),nlx) + nxs
                ird = i
                if (iface(ibc) == 1) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = 2.0d0*Pi*rd(ird-1)*dz_loc_p(ng)
                  
                  !print *,'ptran_conn_bc: ',nc,m,ibc,ng,ird, &
                  !rd(ird-1),areabc(nc)
                  
                  vlbc(nc)   = vlx0
                  if (iphase == 2) vgbc(nc) = vgx0
                else if (iface(ibc) == 2) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = 2.0d0*Pi*rd(ird)*dz_loc_p(ng)
                  
                  !print *,'ptran_conn_bc: ',nc,m,ibc,ng,ird, &
                  !rd(ird-1),areabc(nc)
                  
                  vlbc(nc)   = -vlx0
                  if (iphase == 2) vgbc(nc) = -vgx0
                else if (iface(ibc) == 3) then
                  distbc(nc) = 0.5d0*dz_loc_p(ng)
                  areabc(nc) = Pi*(rd(ird) + rd(ird-1))* &
                               (rd(ird) - rd(ird-1))
                  
                  !print *,'ptran_conn_bc: ',nc,m,ibc,ng,ird, &
                  !rd(ird-1),areabc(nc)
                  
                  vlbc(nc)   = vlz0
                  if (iphase == 2) vgbc(nc) = vgz0
                else if (iface(ibc) == 4) then
                  distbc(nc) = 0.5d0*dz_loc_p(ng)
                  areabc(nc) = Pi*(rd(ird) + rd(ird-1))* &
                               (rd(ird) - rd(ird-1))
                  
                  !print *,'ptran_conn_bc: ',nc,m,ibc,ng,ird, &
                  !rd(ird-1),areabc(nc)
                  
                  vlbc(nc)   = -vlz0
                  if (iphase == 2) vgbc(nc) = -vgz0
                endif
                
                case (3) ! spherical

!               ird = mod(mod((m),nlxy),nlx) + nxs
                ird = i
                if (iface(ibc) == 1) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = 0.d0 ! 4.d0*Pi*rd(0)**2
                  vlbc(nc)   = vlx0
                  if (iphase == 2) vgbc(nc) = vgx0
                else if (iface(ibc) == 2) then
                  distbc(nc) = 0.5d0*dx_loc_p(ng)
                  areabc(nc) = 4.d0*Pi*rd(ird)**2
                  vlbc(nc)   = -vlx0
                  if (iphase == 2) vgbc(nc) = -vgx0
                endif
                
                end select
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    call VecRestoreArrayF90(dx_loc, dx_loc_p, ierr)
    call VecRestoreArrayF90(dy_loc, dy_loc_p, ierr)
    call VecRestoreArrayF90(dz_loc, dz_loc_p, ierr)
    
    if (nconnbc .ne. nc) then
      write(*,*) 'error in computing boundary connections: ', &
      'rank = ',myrank,' nconnbc = ',nconnbc,' nc = ',nc
      stop
    endif

#if 0
    write(*,*) 'ptran_conn: proc',myrank,': nlx=',nlx,' nly=',nly, &
                         ' nlz=',nlz,' ngx=',ngx,' ngy=',ngy,' ngz=',ngz
    write(*,*) 'ptran_conn: proc',myrank,': ','nxs=',nxs, &
    ' nlys=',nlys,'nlzs=',nlzs,' ngxs=',ngxs,' ngys=',ngys,' ngzs=',ngzs
    write(*,*) 'ptran_conn: proc',myrank,': ',nxe,nye,nze,ngxe,ngye, &
                                              ngze
    write(*,*) 'ptran_conn: proc',myrank,': ',istart,iend,jstart,jend, &
                                              kstart,kend,nconn
    write(*,*) 'ptran_conn: proc',myrank,': ',nconnbc,nconnbc_w, &
    nconnbc_e,nconnbc_s,nconnbc_n,nconnbc_t
#endif

    ibug = 0
    if (ibug.eq.1) then
      write(*,*) 'proc nconn  m1 m2 d1 d2 area',nconn
      do nc = 1,nconn
        write(*,*) myrank,nc,nd1(nc),nd2(nc),dist1(nc),dist2(nc),area(nc)
      enddo
      write(*,*) 'proc nconn  mbc dbc areabc',nconnbc
      do nc = 1,nconnbc
        write(*,*) myrank,nc,mblkbc(nc),distbc(nc),areabc(nc)
      enddo
    endif
    
    if (iread_vel == 1) then
      if (myrank == 0) write(*,*) '--> read velocity field'
      call ptran_read_vel (da_1dof)
    endif

  end subroutine ptran_conn
  
!=================================================================

  subroutine ptran_read_vel(da_1dof)

  use ptran_global_module
  use trdynmem_module
  use fileio_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA :: da_1dof
  
  Vec :: N_vec, G_vec, L_vec

  integer :: ierr,i,j,k,fid,m,n,nc
  
  integer, allocatable :: indices(:)

  real*8 :: dum,vel
  
  real*8, pointer :: v_p(:)
  real*8, allocatable :: vel_tmp(:)

  character(len=strlen) :: string2
  character(len=strlen) :: filename

  ierr = 0
  fid = 67
  
!---read FLOTRAN velocity fields (in units of m/y)

    call DACreateNaturalVector(da_1dof,N_vec,ierr)
    call DACreateGlobalVector(da_1dof,G_vec,ierr)
    call DACreateLocalVector(da_1dof,L_vec,ierr)
        
    allocate(vel_tmp(nmax),indices(nmax))

!---x-direction

  if (nx > 1) then
    if (myrank==0) then
    
      filename = flowvx
      open(unit=fid,file=filename,action='read',status='old',iostat=ierr)
      if (ierr /= 0) then
        print *, 'File: ', adjustl(filename), ' not found ', myrank
        stop
      endif

      read(fid,*) string2
      read(fid,*) string2
      read(fid,*) string2

      n = 0
      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            n = n+1
            read(fid,*) dum,dum,dum,vel
            vel_tmp(n) = vel*uyrsec
            indices(n) = n-1
          enddo
          n = n+1
          vel_tmp(n) = -999.d0
          indices(n) = n-1
        enddo
      enddo
      close(fid)
      call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
    endif

    call VecAssemblyBegin(N_vec,ierr)
    call VecAssemblyEnd(N_vec,ierr)

    call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call VecGetArrayF90(L_vec,v_p,ierr)
    
!---set x-velocity for local-ghosted nodes (skip corners and boundaries)
    do nc = 1, nconnx
      m = nd1(nc) ! local/ghosted
      n = nG2L(m) ! local/ghosted -> local/non-ghosted
!     if (n > 0) vl(nc) = v_p(n)
      vl(nc) = v_p(m)
!     print *,'ptranconn-x: ',myrank,nc,m,n,nconnx,nconny,v_p(m)*yrsec
    enddo
    call VecRestoreArrayF90(L_vec,v_p,ierr)
  endif

  if (ny > 1) then
    if (myrank==0 ) then
      filename = flowvy
      open(unit=fid,file=filename,action='read',status='old',iostat=ierr)
      if (ierr /= 0) then
        print *, 'File: ', adjustl(filename), ' not found ', myrank
        stop
      endif
        
      read(fid,*) string2
      read(fid,*) string2
      read(fid,*) string2

      do k=1,nz
        do i=1,nx
          do j=1,ny-1
            n = i+(j-1)*nx+(k-1)*nxy
            read(fid,*) dum,dum,dum,vel
            vel_tmp(n) = vel*uyrsec
            indices(n) = n-1
          enddo
          n = i+(ny-1)*nx+(k-1)*nxy
          vel_tmp(n) = -999.d0
          indices(n) = n-1
        enddo
      enddo
      close(fid)
      call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
    endif

    call VecAssemblyBegin(N_vec,ierr)
    call VecAssemblyEnd(N_vec,ierr)

    call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call VecGetArrayF90(L_vec,v_p,ierr)
    
    do nc = nconnx+1, nconny
      m = nd1(nc)
      n = nG2L(m)
!     if (n > 0) vl(nc) = v_p(m)
      vl(nc) = v_p(m)
!     print *,'ptranconn-y: ',myrank,nc,m,n,nconnx,nconny,v_p(m)*yrsec
    enddo
    call VecRestoreArrayF90(L_vec,v_p,ierr)
  endif
  
!---z-direction

  if (nz > 1) then
    if (myrank == 0) then
      filename = flowvz
      open(unit=fid,file=filename,action='read',status='old',iostat=ierr)
      if (ierr /= 0) then
        print *, 'File: ', adjustl(filename), ' not found ', myrank
      endif
        
      read(fid,*) string2
      read(fid,*) string2
      read(fid,*) string2

        do j=1,ny
          do i=1,nx
            do k=1,nz-1
              n = i+(j-1)*nx+(k-1)*nxy
              read(fid,*) dum,dum,dum,vel
              vel_tmp(n) = vel*uyrsec
              indices(n) = n-1
            enddo
            n = i+(j-1)*nx+(nz-1)*nxy
            vel_tmp(n) = -999.d0
            indices(n) = n-1
          enddo
        enddo
        close(fid)
        call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
      endif

      call VecAssemblyBegin(N_vec,ierr)
      call VecAssemblyEnd(N_vec,ierr)

      call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
      call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
      call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
      call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
      call VecGetArrayF90(L_vec,v_p,ierr)
    
      do nc = nconny+1, nconn
        m = nd1(nc)
        n = nG2L(m)
!       if (n > 0) vl(nc) = v_p(n)
        vl(nc) = v_p(m)
      enddo
      call VecRestoreArrayF90(L_vec,v_p,ierr)
    
  endif

!===============================================

#if 0
!---x-direction

  if (nx > 1) then
    if (myrank==0) then
    
      filename = flowvx
      open(unit=fid,file=filename,action='read',iostat=ierr)
      if (ierr /= 0) then
        print *, 'File: ', adjustl(filename), ' not found ', myrank
        stop
      endif

      read(fid,*) string2
      read(fid,*) string2
      read(fid,*) string2

      n = 0
      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            n = n+1
            read(fid,*) dum,dum,dum,vel
            vel_tmp(n) = vel*uyrsec
            indices(n) = n-1
          enddo
          n = n+1
          vel_tmp(n) = -999.d0
          indices(n) = n-1
        enddo
      enddo
      close(fid)
      call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
    endif

    call VecAssemblyBegin(N_vec,ierr)
    call VecAssemblyEnd(N_vec,ierr)

    call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    
    call VecGetArrayF90(L_vec,v_p,ierr)
    
!---set x-velocity for local-ghosted nodes (skip corners and boundaries)
    nc = 0
    do k = kstart,kend
      do j = jstart,jend
        do i = 1,ngx-1
          ng = i+j*ngx+k*ngxy
          if (v_p(ng) > -999.d0) then
            nc = nc + 1
            vl(nc) = v_p(ng)
          endif
          print *,'ptranconn-x: ',myrank,nc,i,j,k,ng,ngx,ngy,ngxy,v_p(ng)*yrsec
        enddo
      enddo
    enddo
    call VecRestoreArrayF90(L_vec,v_p,ierr)
  endif
    
!---y-direction

  if (ny > 1) then
    if (myrank==0 ) then
      filename = flowvy
      open(unit=fid,file=filename,action='read',iostat=ierr)
      if (ierr /= 0) then
        print *, 'File: ', adjustl(filename), ' not found ', myrank
        stop
      endif
        
      read(fid,*) string2
      read(fid,*) string2
      read(fid,*) string2

      n = 0
!     do k=1,nz
!       do j=1,ny-1
!         do i=1,nx
      do k=1,nz
        do i=1,nx
          do j=1,ny-1
!           n = n + 1
            n = i+(j-1)*nx+(k-1)*nxy
            read(fid,*) dum,dum,dum,vel
            vel_tmp(n) = vel*uyrsec
            indices(n) = n-1
          enddo
!       enddo
!       do i=1,nx
!         n = n + 1
          n = i+(ny-1)*nx+(k-1)*nxy
          vel_tmp(n) = -999.d0
          indices(n) = n-1
!       enddo
        enddo
      enddo
      close(fid)
      call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
    endif

    call VecAssemblyBegin(N_vec,ierr)
    call VecAssemblyEnd(N_vec,ierr)

    call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
    call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
    call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)

    call VecGetArrayF90(L_vec,v_p,ierr)
    
    nc = 0
    do k = kstart,kend
      do i = istart,iend
        do j = 1,ngy-1
!         n = i+1+(j-1)*ngx+k*ngxy
          ng = j+i*ngy+k*ngxy
          if (v_p(ng) > -999.d0) then
            nc = nc + 1
            vl(nc+nconnx) = v_p(ng)
          endif
          print *,'ptranconn-y: ',myrank,nc+nconnx,i,j,k,ng,ngx,ngy,ngxy,v_p(ng)*yrsec
        enddo
      enddo
    enddo
    call VecRestoreArrayF90(L_vec,v_p,ierr)
  endif
  
!---z-direction

    if (nz > 1) then
      if (myrank == 0) then
        filename = flowvz
        open(unit=fid,file=filename,action='read',iostat=ierr)
        if (ierr /= 0) then
          print *, 'File: ', adjustl(filename), ' not found ', myrank
        endif
        
        read(fid,*) string2
        read(fid,*) string2
        read(fid,*) string2

        n = 0
!       do k=1,nz-1
!         do j=1,ny
!           do i=1,nx
        do j=1,ny
          do i=1,nx
            do k=1,nz-1
              n = n + 1
              read(fid,*) dum,dum,dum,vel
              vel_tmp(n) = vel*uyrsec
              indices(n) = n-1
            enddo
!         enddo
!       enddo
!       do j=1,ny
!         do i=1,nx
            n = n + 1
            vel_tmp(n) = -999.d0
            indices(n) = n-1
!         enddo
!       enddo
          enddo
        enddo
        close(fid)
        call VecSetValues(N_vec,nmax,indices,vel_tmp,INSERT_VALUES,ierr)
      endif

      call VecAssemblyBegin(N_vec,ierr)
      call VecAssemblyEnd(N_vec,ierr)

      call DANaturalToGlobalBegin(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
      call DANaturalToGlobalEnd(da_1dof,N_vec,INSERT_VALUES,G_vec,ierr)
      call DAGlobalToLocalBegin(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)
      call DAGlobalToLocalEnd(da_1dof,G_vec,INSERT_VALUES,L_vec,ierr)

      call VecGetArrayF90(L_vec,v_p,ierr)
    
      nc = 0
      do j = jstart,jend
        do i = istart,iend
          do k = 1,ngz-1
            n = i+1+j*ngx+(k-1)*ngxy
            if (v_p(n) > -999.d0) then
              nc = nc + 1
              vl(nc+nconnx+nconny) = v_p(n)
            endif
          enddo
        enddo
      enddo
      call VecRestoreArrayF90(L_vec,v_p,ierr)
    
    endif
#endif

    call VecDestroy(N_vec,ierr)
    call VecDestroy(G_vec,ierr)
    call VecDestroy(L_vec,ierr)
    deallocate(vel_tmp,indices)
    
!   do nc = 1, nconn
!     m = nd1(nc)
!     n = nG2L(m)
!     write(*,*) 'ptranconn-v: ',myrank,nc,m,n,vl(nc)*yrsec
!   enddo

  end subroutine ptran_read_vel

end module ptran_conn_module
