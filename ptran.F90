!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:39:55  lichtner
! Added source term time stepping check.
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

  program ptran

  use ptran_global_module
  use trdynmem_module
  use ptran_read_module
  use ptran_conn_module
  use ptran_init_module
  use ptran_solv_module
  use ptran_dt_module
  use ptran_update_module
  use ptran_out_module
  use ptran_destroy_module

!for running on PC
#ifdef MPICHNT
  use dfwin
#endif
   
  implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#ifdef USE_PETSC216
    ! petsc-2.1.6
#include "include/finclude/petscsles.h"
#endif
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA    :: da, da_mat, da_1dof, da_kin

#ifdef USE_PETSC216
  SLES  :: sles
#endif

#ifdef USE_PETSC221
  integer :: sles
#endif

  PetscLogDouble :: timex(4), timex_wall(4)
  KSP   ::  ksp

!for running on PC
#ifdef MPICHNT
  integer :: prochandle, affmask, retcode  
#endif

  integer :: ierr,i,its,ix,j,jy,kplt=0,kstep,kz,n
  
  real*8 :: tsrc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of Program: PFLOTRAN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ierr = 0
  
! nullify(b_p,phik_p)

! Initialize Startup Time
  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

! Initialize PETSc
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

! Initialize MPI
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)

!for running on PC
#ifdef MPICHNT
  prochandle = GetCurrentProcess()
  if (commsize == 1) then
    affmask = 1
    retcode = SetProcessAffinityMask(prochandle,affmask)
  endif
#endif

  if (myrank==0) then
    print *, '# processors = ', commsize
    print *, ' '
  endif
  
  npx = PETSC_DECIDE
  npy = PETSC_DECIDE
  npz = PETSC_DECIDE
  
  call ptran_read
                   
  call ptran_init (da,da_mat,da_1dof,da_kin,sles,ksp)

  call ptran_chem(da,da_1dof,da_kin,sles)
  
  call ptran_conn (da_1dof)

  if (myrank==0) &
  write(*,*) '--> initialization complete - begin time stepping'
  
  kplt    = 1
  nstep   = 0
  newtcum = 0
  icutcum = 0
  iflgcut = 0
  isrc1 = 2
  t = 0.d0
  do kstep = 1, kmax

    nstep = nstep + 1
    
    t = t + dt

    if (t > tplot(kplt)) then
      t  = t - dt
      dt = tplot(kplt) - t
      t  = tplot(kplt)
    endif
  
!---check source/sink times
    if (nblksrc > 0) then
      n = 1
      tsrc = timesrc(isrc1,n)
      if (t == tsrc) then
        isrc1 = isrc1 + 1
      else if (t > tsrc) then
        t = t - dt
        dt = tsrc - t
        t = tsrc
        isrc1 = isrc1 + 1
!       goto 999
      endif
!999  continue
!     print *,'ptran-src: ',isrc1,tsrc,t,dt
    endif
    
    if (dt <= 0.d0) then
      if (myrank == 0) &
      print *,'ptran-zero dt stopping: ',kstep,t,dt,kplt, &
      tplot(kplt),tsrc
      stop
    endif

!---solve reactive transport equations over a time step
    call ptran_solv (its,da,da_mat,da_1dof,da_kin,sles,ksp)

    newtcum = newtcum + newton
    icutcum = icutcum + icut

!---update field variables
    call ptran_update (da,da_1dof,da_kin)

    if (myrank==0) then
    ! need to check if VecMax is C or FORTRAN numbered!
      imax = imax + 1
!     n = (imax-1)/ncomp+1
      n = imax/ncomp+1
      kz = (n-1)/nxy+1
      jy = (n-1-(kz-1)*nxy)/nx+1
      ix = n-(jy-1)*nx-(kz-1)*nxy
      j = imax-(n-1)*ncomp

      write(*,'(" step   time       dt         dtmax [",a1,"]      cuts", &
  &   "    nwt  iter   R2norm Species")') tunit
      write(*,'(i4,1x,1p3e11.4,1x,i3," [",i6,"]", &
  &   i4," [",i6,"]",1pe12.4,1x,a8)') &
      kstep,t/tconv,dt/tconv,dtmax/tconv, &
      icut,icutcum,newton,newtcum,r2norm,nam(j)

      write(*,'("  dcmax = ",1pe12.4," dc/dt = ",1pe12.4, &
  &   " @ node ",i6," (",i3", ",i3,", ",i3,")")') dcmax0,dcdt,n,ix,jy,kz
      write(*,'("  gmres: ",20i3,(/,10x,20i3))') (itpetsc(i),i=1,newton)
      do i = 1, newton
        itpetscum = itpetscum + itpetsc(i)
      enddo
      write(*,*)
    endif

    ! RTM: Adjust the timestep.
    call ptran_dt (nstep,newton,t,dt,dtmax,tfac,iflgcut,iaccel,myrank)
    
    call ptran_psai_out (kplt,da,da_1dof,da_kin)
    if (t >= tplot(kplt)) kplt = kplt + 1
    if (kplt > kplot) exit
  enddo

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  call PetscGetTime(timex_wall(2), ierr)
        
      if (myrank == 0) then
        write(*,'("steps: ",13x,i6,/, &
    &   "total newton iters: ",i6,/, &
    &   "total PETSc gmres:  ",i6,/, &
    &   "total cuts: ",11x,i3)') nstep,newtcum,itpetscum,icutcum

        write(*,'("CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(*,'("Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0

        write(iunit2,'(/,"steps: ",13x,i6,/, &
    &   "total newton iters: ",i6,/, &
    &   "total PETSc gmres:  ",i6,/, &
    &   "total cuts: ",11x,i3)') nstep,newtcum,itpetscum,icutcum

        write(iunit2,'("CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(IUNIT2,'("Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0
        
        close (iunit2)
        close (ptran_iunit4)
      endif
! -----------------------------------------

  call ptran_destroy (da,da_mat,da_1dof,da_kin,sles)

  call PetscFinalize(ierr)

end program ptran
