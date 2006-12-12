!=======================================================================
! PFLOTRAN v1.0 LA-CC-06-093
!=======================================================================

! The software titled "PFLOTRAN v 1.0" has been assigned LA-CC-06-093. 
! The software is unclassified and does not contain Unclassified 
! Controlled Nuclear Information (UCNI).

! The software is under review to be released as open source, which 
! requires review and approval by the appropriate DOE Program Office. 
! If DOE declines to release this software as publicly available open 
! source, it would be subject to export control under Department of 
! Commerce regulations, and classified as ECCN 
! (Export Control Classification Number) EAR99. 
 
! If released as open source, the software will be publicly available and 
! not subject to export control. Until DOE approval is obtained, the 
! software should be treated as export controlled. DOE reserves the 
! right to release or deny release of the software as open source.

! Send all bug reports/questions/comments to:
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-6, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM
!=======================================================================
  program pflotran
  use pflow_gridtype_module
  use pflow_grid_module
  use pflow_read_gridsize_module
  use pflow_output_module
  use span_wagner_module
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
  use rock_reac_module
  use pflotran_couple_module
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscpc.h"
!#include "include/finclude/petscsles.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  integer :: ierr
  type(pflowGrid) :: grid
! integer :: nx, ny, nz
! integer :: npx, npy, npz
  integer :: nphase, ndof, icouple, itable
  integer :: nspec,npricomp
  integer :: kplt, iplot, igeom, iflgcut_pflow, its_pflow, ntstep
  integer :: myid
  integer :: steps
  
  real*8 :: tflow1, tflow2, sumdt
! real*8 ::  tmp1, tmp2, drtotmax, tmp3,tmp4

!==========trans==========
  DA    :: da, da_mat, da_1dof, da_kin

#ifdef USE_PETSC216
  SLES  :: sles
#endif

#ifdef USE_PETSC221
  integer :: sles
#endif

  KSP   :: ksp
  integer :: i,its,ix,j,jy,kstep,kz,n,ihalcnt
!==========trans==========

! Initialize Startup Time
  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

! Initialize PETSc
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

! Initialize MPI
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, myid, ierr)

!======================= initialize pflow =======================
  call pflow_read_gridsize("pflow.in", igeom, nx, ny, nz, npx, npy, npz, &
  nphase, nspec, npricomp, ndof, idcdm, itable, ierr)
  
  icouple = 1

  if(ierr /= 0) then
    if(myid == 0) then
      print *, "Error reading gridsize--exiting."
    endif 
    call PetscFinalize(ierr)
    stop
  endif

! set up structure constructor
! npx = PETSC_DECIDE; npy = PETSC_DECIDE; npz = PETSC_DECIDE
  grid = pflowGrid_new(igeom, nx, ny, nz, npx, npy, npz, nphase, nspec, npricomp, &
         ndof, icouple, idcdm, itable)
  if (nphase>=2) call initialize_span_wagner(itable,grid%myrank)
  
  call pflowGrid_setup(grid, "pflow.in")

  kplt = 0
  iplot = 1
  
  call pflow_output(grid,kplt,iplot)

  if(myid == 0) print *, ""
  if(myid == 0) print *, ""

!======================= initialize ptran =======================

  npx = PETSC_DECIDE
  npy = PETSC_DECIDE
  npz = PETSC_DECIDE

  call ptran_read
  print *, 'finished  read'
  call ptran_init(da,da_mat,da_1dof,da_kin,sles,ksp)
  print *, 'finished  init'
  call pflowGrid_ptran_init(grid, ierr)
  print *, 'finished  flow init'
  call ptran_chem(da,da_1dof,da_kin,sles)
  print *, 'finished  chem'
  if (grid%ihydrostatic == 1) call ptran_bc_reassign(grid)
  print *, 'finished  rebc'
!  call ptran_conn(da_1dof)

  if (myrank==0) &
  write(*,*) '--> initialization complete - begin time stepping'

!======================= begin time stepping =======================
  kplt = 1
  iplot = 0
  grid%isrc1 = 2
  
  nstep   = 0
  newtcum = 0
  icutcum = 0
  iflgcut_pflow = 0
  iflgcut = 0
  t = 0.d0
  tflow2 = 0.d0
  rtoth2o = 0.D0
  rtotco2 = 0.D0
  grid%rtot = 0.D0 

! flow time step
  do steps = 1, grid%stepmax
    tflow1 = tflow2

!   pass reaction rates for minerals containing CO2 and H2O
    grid%rtot(:,1)= rtoth2o(:)
    grid%rtot(:,2)= rtotco2(:)
!   print *,'Pflow: rtot=', grid%rtot(1:128,1), grid%rtot(1:128,2) 
    call VecCopy(porosity,grid%porosity,ierr)

    call pflowGrid_step(grid,ntstep,kplt,iplot,iflgcut_pflow,ihalcnt,its_pflow)
    
!   print *,'pflotran: grid%vvlbc ', grid%vvlbc
!   print *,'pflotran: grid%vvlbc ',grid%vvgbc
    
!   compute velocity at time tflow2
!   call pflowgrid_setvel(grid, grid%vvl_loc, grid%vvlbc, ibconn, ibndtyp)

    tflow2 = grid%t
    t = tflow1
    if (grid%isync == 1) dt = grid%dt

!   trans time step
    rtoth2o = 0.D0
    rtotco2 = 0.D0
    sumdt = 0.d0
    do kstep = 1, kmax

      nstep = nstep + 1

      t = t + dt

      if (1.2d0 * t > tflow2) then
        t  = t - dt
        dt = tflow2 - t
        t  = tflow2
      endif
      
!     interpolate field variables
      call pflotranGrid_interp(grid, tflow1, tflow2)
 
!     solve reactive transport equations over a time step
      call ptran_solv(its,da,da_mat,da_1dof,da_kin,sles,ksp)

      newtcum = newtcum + newton
      icutcum = icutcum + icut

!     update field variables
      call ptran_update(da,da_1dof,da_kin)
      if(ipor>0) call Rock_reac(grid)

!     time average rtot
      rtoth2o = rtoth2o + rrtoth2o*dt
      rtotco2 = rtotco2 + rrtotco2*dt
		
!     drtotmax=0.D0; tmp3=0.D0; tmp4=0.D0	
!     do n=1,grid%nlmax
!       tmp1=dabs(rtoth2o(n)-rrtoth2o(n))
!       tmp2=dabs(rtotco2(n)-rrtotco2(n))
!       if(tmp3<dabs(rrtoth2o(n))) tmp3=dabs(rrtoth2o(n))
!       if(tmp4<dabs(rrtotco2(n))) tmp4=dabs(rrtotco2(n))
!       if(drtotmax<tmp1) drtotmax=tmp1
!       if(drtotmax<tmp2) drtotmax=tmp2
!     enddo
!     PRINT *, 'pflotran/rtot :: ',tmp3,tmp4,drtotmax,rrtoth2o(3)

      if (myrank==0) then
        sumdt = sumdt + dt
!       print *,'pflotran: ',tflow2,tflow1,tflow2-tflow1,sumdt,dt
        
        ! need to check if VecMax is C or FORTRAN numbered!
        imax = imax + 1
!       n = (imax-1)/ncomp+1
        n = imax/ncomp+1
        kz = (n-1)/nxy+1
        jy = (n-1-(kz-1)*nxy)/nx+1
        ix = n-(jy-1)*nx-(kz-1)*nxy
        j = imax-(n-1)*ncomp
        
        if (j == 0) j = 1

        write(*,'(/," TRANS     time       dt         target [",a1,"]     cuts", &
  &     "      nwt iter     R2norm Species")') tunit
        write(*,'(i2,"/",i4,1x,1p3e11.4,1x,i3," [",i6,"]", &
  &     i4," [",i6,"]",1pe12.4,1x,a8)') &
        kstep,nstep,t/tconv,dt/tconv,tflow2/tconv, &
        icut,icutcum,newton,newtcum,r2norm,nam(j)

        write(*,'("  dcmax = ",1pe12.4," dc/dt = ",1pe12.4, &
  &     " @ node ",i6," (",i3", ",i3,", ",i3,")")') dcmax0, &
        dcdt,n,ix,jy,kz
        write(*,'("  gmres: ",20i3,(/,10x,20i3))') &
        (itpetsc(i),i=1,newton)
        do i = 1, newton
          itpetscum = itpetscum + itpetsc(i)
        enddo
        write(*,'(80("="))')
      endif

      ! Adjust ptran timestep
      if (grid%isync == 0) &
      call ptran_dt (nstep,newton,t,dt,dtmax,tfac,iflgcut,iaccel,myrank)

      call ptran_psi_out(kplt,da,da_1dof,da_kin)
      
      if (t >= tflow2) exit
    enddo

!   update pflow field variables
    call pflowgrid_update(grid)

    rtoth2o = rtoth2o/(tflow2-tflow1)
    rtotco2 = rtotco2/(tflow2-tflow1)
!   print *,'Ptran: rtot=', rtoth2o(1:128),rtotco2(1:128)

    call pflow_output(grid,kplt,iplot)
    if (iflgcut_pflow == 0) call pflowgrid_update_dt(grid,its_pflow)
    
!   print *,'pflotran: ',steps,iflgcut,grid%t,grid%dt,its
!   print *,steps,grid%kplot,grid%t,grid%tplot(grid%kplot)
    if (grid%t >= grid%tplot(grid%kplot)) exit  
  enddo

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  call PetscGetTime(timex_wall(2), ierr)
        
  if (myrank == 0) then

!   write(*,'("steps: ",13x,i6,/, &
!   &   "total newton iters: ",i6,/, &
!   &   "total PETSc gmres:  ",i6,/, &
!   &   "total cuts: ",11x,i3)') steps,newtcum,itpetscum,icutcum
        
    write(*,'(/,"PFLOW steps:        ",i6)') steps
    write(IUNIT2,'(/,"PFLOW steps:        ",i6)') steps

    write(*,'("+++++++++++++++++++++++++++")')
    write(IUNIT2,'("+++++++++++++++++++++++++++")')
        
    write(*,'( &
    &   "PTRAN steps:        ",i6,/, &
    &   "total newton iters: ",i6,/, &
    &   "total PETSc gmres:  ",i6,/, &
    &   "total cuts: ",11x,i3)') nstep,newtcum,itpetscum,icutcum

    write(iunit2,'(/,"PTRAN steps: ",13x,i6,/, &
    &   "total newton iters: ",i6,/, &
    &   "total PETSc gmres:  ",i6,/, &
    &   "total cuts: ",11x,i3)') nstep,newtcum,itpetscum,icutcum

    write(*,'("+++++++++++++++++++++++++++")')
    write(IUNIT2,'("+++++++++++++++++++++++++++")')

    write(*,'("Total CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

    write(*,'("Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0

    write(iunit2,'("Total CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

    write(iunit2,'("Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0
        
  endif

  call pflowgrid_destroy(grid)

  call ptran_destroy(da,da_mat,da_1dof,da_kin,sles)

  if (myrank == 0) then
    close(IUNIT2)
    if (grid%ibrkcrv > 0) close(IUNIT4)
    if (ibrkcrv > 0) close(ptran_IUNIT4)
  endif

  call PetscFinalize(ierr)
  
  end program pflotran
