
module pflow_output_module

  public
  
 integer,private :: size_var_use, size_var_node

contains

subroutine pflow_output(grid,kplt,iplot)
  
  use pflow_gridtype_module
  use TTPHASE_module
  use PetscRelWrappers  ! For petsc-release compatibility.
  
  implicit none

#include "include/finclude/petsc.h"
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(inout) :: iplot, kplt
   
   
  Vec :: p_nat, t_nat, c_nat, phis_nat, por_nat, vl_nat, s_nat, x_nat, perm_nat,&
         iphase_nat
  Vec :: c_all, vl_all, p_all, t_all, s_all, x_all, phis_all, iphase_all
  VecScatter :: scat_1dof, scat_3npdof, scat_nph

  real*8, pointer :: p_p(:), t_p(:), c_p(:), phis_p(:), vl_p(:), s_p(:)
  real*8, pointer :: x_p(:), iphase_p(:)
  real*8, pointer :: fldflx(:), fldvol(:)
  
  real*8 :: vavel(3),vaveg(3)
  integer :: nxm1,nym1,nzm1
  
  integer :: ierr
  integer :: i,ip,j,jn,k,ibrk,n,nc,nn
  
  real*8 :: tyr, sum1, sum2, sum1v, sum2v, area, vol, vel, vf
  integer, save :: icall, icall_brk, icall_plt
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab

  PetscTruth :: use_soutput
  
  data icall/1/, icall_brk/1/, icall_plt/1/
  
   
  !save icall, icall_brk, icall_plt
  
! ibrkcrv = 0, no time-history
!         = 1, time-history output
  
! iplot  = -1, Glenn's format
! iplot  = 0, only time-history
!        = 1, spatial plot

! iprint = -1, none
!        =  0, p, t, s, c
!        =  1, vl
!        =  2, porosity, (permeability)
  !return

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_soutput", &
                           use_soutput, ierr)
  if (use_soutput == PETSC_TRUE) then
    call pflow_soutput(grid, kplt, iplot)
    return
  endif
  
  if (iplot == 1 .and. grid%iprint == -2) then
    call geh_io(grid,kplt)
!    kplt = kplt + 1
!    iplot = 0
!    return
  endif
  
  if ((grid%ibrkcrv == 0 .and. iplot == 0) .or. grid%iprint == -1) then
    if (grid%iprint==-1 .and. iplot==1) then
      kplt = kplt + 1
      iplot = 0
    endif
    return
  endif
  
!  if(grid%nphase>1) call pflow_2phase_massbal(grid)
      
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1 .and. icall == 1) then
    icall = 0
    if (grid%myrank == 0) then
      write(fname,'(a9,a4)') 'pflow_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' iplot= ',iplot
      open(unit=IUNIT4,file=fname,action="write")
      write(IUNIT4,'("%#t            dt",100i12)') (i,i=1,grid%ibrkcrv), &
      (i,i=1,grid%ibrkcrv)
    endif
    if (iplot == 0) return
  endif
  
  if (grid%ndof == 1 .and. iplot == 0) return ! no time-history plot for 1 dof

! Create Natural Vec for output: use VecDuplicate here?
  call DACreateNaturalVector(grid%da_1_dof, c_nat, ierr)
  call VecDuplicate(c_nat, phis_nat, ierr)
  call VecDuplicate(c_nat, t_nat,    ierr)
  call VecDuplicate(c_nat, por_nat,  ierr)
  call VecDuplicate(c_nat, perm_nat,  ierr)
  call VecDuplicate(c_nat, iphase_nat, ierr)
  call DACreateNaturalVector(grid%da_nphase_dof, p_nat, ierr)
  call VecDuplicate(p_nat, s_nat, ierr)
  call VecDuplicate(p_nat, x_nat, ierr)
  call DACreateNaturalVector(grid%da_3np_dof, vl_nat, ierr)
  
  ! time-history plots
  
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1) then

!   concentration
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat, &
                                ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat,ierr)

!   velocity fields
    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat, &
                                ierr)
    call DAGlobalToNaturalEnd(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must be
! followed by a call to VecDestroy(x_all)


    call VecScatterCreateToAll(c_nat, scat_1dof, c_all, ierr)
    call VecScatterBegin(scat_1dof, c_nat, c_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_1dof, c_nat, c_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)

    call VecScatterCreateToAll(vl_nat, scat_3npdof, vl_all, ierr)
    call VecScatterBegin(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
       
    call VecScatterDestroy(scat_1dof, ierr)
    call VecScatterDestroy(scat_3npdof, ierr)

    if (grid%myrank == 0) then
      call VecGetArrayF90(c_all, c_p, ierr)
      call VecGetArrayF90(vl_all, vl_p, ierr)
      
      allocate(fldflx(grid%ibrkcrv))
      allocate(fldvol(grid%ibrkcrv))
      
      fldflx = 0.d0
      fldvol = 0.d0
      do ibrk = 1, grid%ibrkcrv
        sum1 = 0.d0
        sum2 = 0.d0
        sum1v = 0.d0
        sum2v = 0.d0
        do k = grid%k1brk(ibrk),grid%k2brk(ibrk)
          do j = grid%j1brk(ibrk),grid%j2brk(ibrk)
            do i = grid%i1brk(ibrk),grid%i2brk(ibrk)
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
              vol = grid%dx0(i) * grid%dy0(j) * grid%dz0(k)
              if (grid%ibrkface(ibrk) == 1) then
                nn = i-1+(j-1)*grid%nx+(k-1)*grid%nxy
                area = grid%dy0(j) * grid%dz0(k)
              else if (grid%ibrkface(ibrk) == 2) then
                nn = i+(j-2)*grid%nx+(k-1)*grid%nxy
                area = grid%dx0(i)*grid%dz0(k)
              else if (grid%ibrkface(ibrk) == 3) then
                nn = i+(j-1)*grid%nx+(k-2)*grid%nxy
                area = grid%dx0(i)*grid%dy0(j)
              endif
              
              vel = vl_p(grid%ibrkface(ibrk)+3*(nn-1))

              sum1 = sum1 + vel*area*c_p(n)
              sum2 = sum2 + vel*area
              sum1v = sum1v + vol*c_p(n)
              sum2v = sum2v + vol
              
!             print *,'output-brk: ',grid%ibrkcrv,ibrk,i,j,k,n,nn, &
!             grid%ibrkface(ibrk),area, &
!             vol,vel*grid%tconv,c_p(n),sum1,sum2,sum1v,sum2v
            enddo
          enddo
        enddo
        if (sum2 .ne. 0.d0) fldflx(ibrk) = sum1/sum2
        if (sum2v .ne. 0.d0) fldvol(ibrk) = sum1v/sum2v
      enddo

      write(IUNIT4,'(1p100e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv, &
                                   (fldflx(i),i=1,grid%ibrkcrv), &
                                   (fldvol(i),i=1,grid%ibrkcrv)
      
      call VecRestoreArrayF90(c_all, c_p, ierr)
      call VecRestoreArrayF90(vl_all, vl_p, ierr)
      deallocate(fldflx)
      deallocate(fldvol)
    endif
    call VecRestoreArrayF90(c_all, c_p, ierr)
    call VecRestoreArrayF90(vl_all, vl_p, ierr)
    call VecDestroy(c_all, ierr)
    call VecDestroy(vl_all, ierr)

    if (iplot == 0)then
      call VecDestroy( c_nat, ierr)
      call VecDestroy( iphase_nat, ierr)
      call VecDestroy( phis_nat, ierr)
      call VecDestroy( t_nat,    ierr)
      call VecDestroy( por_nat,  ierr)
      call VecDestroy( p_nat, ierr)
      call VecDestroy( s_nat, ierr)
      call VecDestroy( x_nat, ierr)
      call VecDestroy( vl_nat, ierr)
      return
    endif
  endif
  
! plot spatial data (iplot > 0)
  
! presssure field
  call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%pressure,INSERT_VALUES, &
                              p_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%pressure,INSERT_VALUES, &
                            p_nat,ierr)

! temperature field
  call DAGlobalToNaturalBegin(grid%da_1_dof,grid%temp,INSERT_VALUES, &
                              t_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_1_dof,grid%temp,INSERT_VALUES, &
                            t_nat,ierr)

! saturation
  call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%sat,INSERT_VALUES, &
                              s_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%sat,INSERT_VALUES, &
                            s_nat,ierr)
  
! primary variables
  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE) then
    call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%xmol,INSERT_VALUES, &
                                x_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%xmol,INSERT_VALUES, &
                              x_nat,ierr)
  endif  
  
  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%iphas,INSERT_VALUES, &
                                iphase_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%iphas,INSERT_VALUES, &
                              iphase_nat,ierr)
  endif  
  
! concentration: check if already done in time-history plot
  if (grid%ibrkcrv == 0) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%conc,INSERT_VALUES, &
                                c_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%conc,INSERT_VALUES, &
                              c_nat,ierr)
  endif
  
! call VecView(phis,PETSC_VIEWER_STDOUT_WORLD,ierr)

! solid volume fraction
  if (grid%rk > 0.d0) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%phis,INSERT_VALUES, &
                                phis_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%phis,INSERT_VALUES, &
                              phis_nat,ierr)
  endif
  
  icall_plt = 0
#ifdef HAVE_MPITOMPIZERO
  call VecScatterCreateToZero(p_nat, scat_nph,  p_all, ierr)
  if (grid%ibrkcrv == 0) &
    call VecScatterCreateToZero(c_nat,  scat_1dof, c_all, ierr)
!   call VecScatterCreateToZero(vl_nat,   scat_3dof, vl_all, ierr)
  call VecScatterCreateToZero(s_nat, scat_nph,  s_all, ierr)
  call VecScatterCreateToZero(t_nat, scat_1dof, t_all, ierr)
  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE ) &
    call VecScatterCreateToZero(x_nat, scat_nph, x_all, ierr)

  if (rk > 0.d0) &
    call VecScatterCreateToZero(phis_nat, scat_1dof, phis_all, ierr)
!   call VecScatterCreateToZero(por_nat,  scat_1dof, por_all, ierr)
#else
  call VecScatterCreateToAll(p_nat, scat_nph,  p_all, ierr)
  call VecScatterCreateToAll(c_nat, scat_1dof, c_all, ierr)
! call VecScatterCreateToAll(vl_nat, scat_3dof, vl_all, ierr)
! call VecScatterCreateToAll(s_nat, scat_nph,  s_all, ierr)
! call VecScatterCreateToAll(t_nat, scat_1dof, t_all, ierr)

  call VecDuplicate(p_all,s_all,ierr)
  call VecDuplicate(c_all,t_all,ierr)
    
  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE) &
    call VecDuplicate(p_all,x_all,ierr)

  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
      .or. grid%use_flash == PETSC_TRUE .or. grid%use_richards == PETSC_TRUE) &
    call VecDuplicate(c_all,iphase_all,ierr)

  if (grid%rk > 0.d0) then
!   call VecScatterCreateToAll(phis_nat, scat_1dof, phis_all, ierr)
    call VecDuplicate(c_all,phis_all,ierr)
  endif
!   call VecScatterCreateToAll(por_nat,  scat_1dof, por_all, ierr)
#endif

! note that the following routines call VecCreate(x_all) and therefore must be
! followed by a call to VecDestroy(x_all)

#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(p_nat, p_all, ierr);   CHKERRQ(ierr)
! call VecConvertMPIToMPIZero(t_nat, t_all, ierr);   CHKERRQ(ierr)
! call VecConvertMPIToMPIZero(s_nat, s_all, ierr);   CHKERRQ(ierr)

! call VecScatterCreateToZero(p_nat, scat_nph, p_all, ierr)
  call VecScatterBegin(scat_nph, p_nat, p_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  call VecScatterEnd(scat_nph, p_nat, p_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

! call VecScatterCreateToZero(s_nat, scat_nph, s_all, ierr)
  call VecScatterBegin(scat_nph, s_nat, s_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  call VecScatterEnd(scat_nph, s_nat, s_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE
       .or. grid%use_richards == PETSC_TRUE) then
    call VecScatterBegin(scat_nph, x_nat, x_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_nph, x_nat, x_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  endif
  
! call VecScatterCreateToZero(t_nat, scat_1dof, t_all, ierr)
  call VecScatterBegin(scat_1dof, t_nat, t_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, t_nat, t_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

  if (ibrkcrv == 0) then
!   call VecConvertMPIToMPIZero(c_nat, c_all, ierr)

!   call VecScatterCreateToZero(c_nat, scat_1dof, c_all, ierr)
    call VecScatterBegin(scat_1dof, c_nat, c_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_1dof, c_nat, c_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  endif
  
  if (rk > 0.d0) then
!   call VecConvertMPIToMPIZero(phis_nat, phis_all, ierr)

!   call VecScatterCreateToZero(phis_nat, scatter, phis_all, ierr)
    call VecScatterBegin(scat_1dof, phis_nat, phis_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_1dof, phis_nat, phis_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  endif
#else
! call VecConvertMPIToSeqAll(p_nat, p_all, ierr);    CHKERRQ(ierr)
! call VecConvertMPIToSeqAll(t_nat, t_all, ierr);    CHKERRQ(ierr)
! call VecConvertMPIToSeqAll(s_nat, s_all, ierr);    CHKERRQ(ierr)

! call VecScatterCreateToAll(p_nat, scatter, p_all, ierr)
  call VecScatterBegin(scat_nph, p_nat, p_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  call VecScatterEnd(scat_nph, p_nat, p_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

! call VecScatterCreateToAll(s_nat, scatter, s_all, ierr)
  call VecScatterBegin(scat_nph, s_nat, s_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  call VecScatterEnd(scat_nph, s_nat, s_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE) then
    call VecScatterBegin(scat_nph, x_nat, x_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_nph, x_nat, x_all, INSERT_VALUES, SCATTER_FORWARD, &
                       ierr)
  endif
  
  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
      .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE) then
    call VecScatterBegin(scat_1dof, iphase_nat, iphase_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_1dof, iphase_nat, iphase_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)

  endif
! call VecScatterCreateToAll(t_nat, scatter, t_all, ierr)
  call VecScatterBegin(scat_1dof, t_nat, t_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, t_nat, t_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)

!  if (ibrkcrv == 0) then
!   call VecConvertMPIToSeqAll(c_nat, c_all, ierr)

!   call VecScatterCreateToAll(c_nat, scat_1dof, c_all, ierr)
  call VecScatterBegin(scat_1dof, c_nat, c_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, c_nat, c_all, INSERT_VALUES, SCATTER_FORWARD, &
                     ierr)
!  endif

  if (grid%rk > 0.d0) then
!   call VecConvertMPIToSeqAll(phis_nat, phis_all, ierr)

!   call VecScatterCreateToAll(phis_nat, scat_1dof, phis_all, ierr)
    call VecScatterBegin(scat_1dof, phis_nat, phis_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_1dof, phis_nat, phis_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  endif

  call VecScatterDestroy(scat_nph,  ierr)
  call VecScatterDestroy(scat_1dof, ierr)

#endif

  tab = char(9)
  q = '","'
  tyr = grid%t/grid%tconv

  if (grid%myrank == 0) then
    if (grid%flowsteps == 0) then
      call SNESGetTolerances(grid%snes, grid%atol, grid%rtol, grid%stol, &
                             grid%maxit, grid%maxf, ierr)
!       write(*,'("%#atol, rtol, stol= ",1p3e12.4," maxit, maxf= ",2i6)') &
!       grid%atol,grid%rtol,grid%stol,grid%maxit,grid%maxf
  
!     Calculate the x, y, z vectors that give the 
!     physical coordinates of each cell.
!     call pflowGrid_compute_xyz(grid)
    endif
    
    call VecGetArrayF90(p_all, p_p, ierr)
    call VecGetArrayF90(t_all, t_p, ierr)
    call VecGetArrayF90(c_all, c_p, ierr)
    call VecGetArrayF90(s_all, s_p, ierr)
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
        grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
         .or. grid%use_richards == PETSC_TRUE )&
      call VecGetArrayF90(x_all, x_p, ierr)


    if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
        .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)&
      call VecGetArrayF90(iphase_all, iphase_p, ierr)

    if (grid%rk > 0.d0) then
      call VecGetArrayF90(phis_all, phis_p, ierr)
    endif
    
    if (kplt < 10) then
      write(fname,'(a7,i1,a4)') 'pflow00', kplt, '.dat'
    else if (kplt < 100) then
      write(fname,'(a6,i2,a4)') 'pflow0', kplt, '.dat'
    else
      write(fname,'(a5,i3,a4)') 'pflow', kplt, '.dat'
    endif
    write(*,*) '--> write output file: ',fname
    open(unit=IUNIT3,file=fname,action="write")
      
    if (grid%nx+grid%ny+grid%nz .gt. grid%nx*grid%ny*grid%nz) then ! 1D
    
      write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," [",a1, &
        & "]", " dt =",1pe12.4," sec")') &
        grid%flowsteps,grid%t,tyr,grid%tunit,grid%dt
      if (grid%use_2ph == PETSC_TRUE) then
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
          & "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
          & "-xlco2-     -xgco2-        -Vf-     ")')
      elseif( grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
              .or. grid%use_flash == PETSC_TRUE) then 
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   iphase ", &
          & "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
          & "-xlco2-     -xgco2-        -Vf-     ")')
      else
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
          & "  -p[Pa]-      -t[C]-      -sl(g)-      -C[mol/L]-      -Vf- &
          & ")')
      endif

      do n = 1, grid%nmax
        jn = grid%jh2o+(n-1)*grid%nphase
        vf = 0.d0
        if (grid%rk > 0.d0) vf = phis_p(n)
        if (grid%use_2ph == PETSC_TRUE ) then
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn),s_p(jn+1), &
            x_p(jn),x_p(jn+1),vf
        elseif (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
                .or. grid%use_flash == PETSC_TRUE)then !.or. grid%use_richards == PETSC_TRUE)then
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            iphase_p(n),&
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn:jn+grid%nphase-1), &
            x_p(jn:jn+grid%nphase-1),vf
      !    print *,'output', grid%nphase, n, (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn:jn+grid%nphase-1), &
      !      x_p(jn:jn+grid%nphase-1),vf
        else
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o), t_p(n), s_p(jn), c_p(n), vf
        endif
      enddo

    else if (grid%nz == 1) then ! 2D x-y

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
          grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE)then
  
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
      else
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
      endif
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
        & '' , J='',i4)') tyr,grid%nx,grid%ny

      do n = 1, grid%nmax
        jn = grid%jh2o+(n-1)*grid%nphase
        vf = 0.d0
        if (grid%rk > 0.d0) vf = phis_p(n)
        if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
            grid%use_vadose == PETSC_TRUE.or.grid%use_flash == PETSC_TRUE &
             ) then
          write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n), &
             p_p(jn+1), t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
        else
          write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), &
            p_p(jn), t_p(n), s_p(jn), c_p(n), vf
        endif
      enddo
          
    else if (grid%ny == 1) then ! 2D x-z

      if (grid%itecplot == 0) then
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        if( grid%use_2ph == PETSC_TRUE ) then
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
         elseif(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
                .or. grid%use_flash == PETSC_TRUE )then
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
        else
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
        endif
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
          & '' , K='',i4)') tyr,grid%nx,grid%nz
        
        do n = 1, grid%nmax
          jn = grid%jh2o+(n-1)*grid%nphase
          vf = 0.d0
          if (grid%rk > 0.d0) vf = phis_p(n)
          if ( grid%use_2ph == PETSC_TRUE) then
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn+1), t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
          elseif(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
                 .or. grid%use_flash == PETSC_TRUE &
                 ) then
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n),iphase_p(n), &
              p_p(jn+1),t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
          else
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn), t_p(n), s_p(jn), c_p(n), vf
          endif
        enddo

      else

        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
            grid%use_vadose == PETSC_TRUE.or.grid%use_flash == PETSC_TRUE &
             ) then
           write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
             'x',q,'z',q,'phase',q,'pl',q,'pg',q,'T',q,'sl',q,'sg',q,'xl',q, &
             'xg',q,'vf','"'
        else
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
        endif
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
          & '' , J='',i4)') tyr,grid%nx,grid%nz
        
!        do i = 1, grid%nx
!          do k = 1, grid%nz
!            n = i+(k-1)*grid%nxy
        do n=1, grid%nmax
          jn = grid%jh2o+(n-1)*grid%nphase
          vf = 0.d0
          if (grid%rk > 0.d0) vf = phis_p(n)
            
!            print *,'pflow_output: ',n,jn,p_p(jn),p_p(jn+1), t_p(n), &
!                                    s_p(jn),s_p(jn+1), x_p(jn),x_p(jn+1)
            
          if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
              .or. grid%use_vadose == PETSC_TRUE &
              .or. grid%use_flash == PETSC_TRUE &
              ) then
            write(IUNIT3,'(1p100e12.4)') grid%x(n),grid%z(n),iphase_p(n), &
              p_p(jn),p_p(jn+1),t_p(n),s_p(jn),s_p(jn+1),x_p(jn),x_p(jn+1),vf
          else
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn), t_p(n), s_p(jn), c_p(n), vf
          endif
!          enddo
        enddo
      endif

    else ! 3D

    write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
    if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   ) then
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'y',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
    else
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'y',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
    endif
    write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
      & '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
    do n = 1,  grid%nmax
      jn = grid%jh2o+(n-1)*grid%nphase
      vf = 0.d0
      if (grid%rk > 0.d0) vf = phis_p(n)
      if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
          .or. grid%use_flash == PETSC_TRUE &   
          .or. grid%use_vadose == PETSC_TRUE) then
        write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), grid%z(n), &
          iphase_p(n), p_p(jn+1), t_p(n), s_p(jn+1), x_p(jn),x_p(jn+1), vf
      else
        write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), grid%z(n), &
          p_p(jn), t_p(n), s_p(jn), c_p(n), vf
      endif
    enddo
  endif

!   write out initial conditions

  if (grid%write_init == 1) then
    fname = 'pflow_init0.dat'
    write(*,*) '--> write output file: ',fname
    close (IUNIT3)
    open(unit=IUNIT3,file=fname,action="write")
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
        .or. grid%use_flash == PETSC_TRUE & 
        .or. grid%use_vadose == PETSC_TRUE &
        .or. grid%use_richards == PETSC_TRUE) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2", &
          & "       p      ","      T      ","     sl(g)      ", &
          & "       xl       xg   ")')
      else       
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2", &
          & "       p      ","      T      ","     sl(g)      ", &
          & "      C      ")')
      endif
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
            jn = grid%jh2o+(n-1)*grid%nphase
            if (grid%use_2ph == PETSC_TRUE .or. &
                grid%use_mph == PETSC_TRUE .or. &
                grid%use_flash == PETSC_TRUE .or. &
                grid%use_vadose == PETSC_TRUE .or.&
                grid%use_richards == PETSC_TRUE ) then
              write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') i,i,j,j,k,k, &
                p_p(jn), t_p(n), s_p(jn), x_p(jn),x_p(jn+1)
            else
              write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') i,i,j,j,k,k, &
                p_p(jn), t_p(n), s_p(jn), c_p(n)
            endif
          enddo
        enddo
      enddo
      write(IUNIT3,'("/",/)')
    endif
    
    call VecRestoreArrayF90(p_all, p_p, ierr)
    call VecRestoreArrayF90(t_all, t_p, ierr)
    call VecRestoreArrayF90(c_all, c_p, ierr)
    call VecRestoreArrayF90(s_all, s_p, ierr)
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
        .or. grid%use_flash == PETSC_TRUE &
        .or. grid%use_vadose == PETSC_TRUE &
         .or. grid%use_richards == PETSC_TRUE ) then
      call VecRestoreArrayF90(x_all, x_p, ierr)
      call VecRestoreArrayF90(iphase_all, iphase_p, ierr)
    endif
    if (grid%rk > 0.d0) then
      call VecRestoreArrayF90(phis_all, phis_p, ierr)
    endif
  endif
  
  call VecDestroy(p_all, ierr)
  call VecDestroy(t_all, ierr)
  call VecDestroy(c_all, ierr)
  call VecDestroy(s_all, ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
      .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE ) &
    call VecDestroy(x_all, ierr)
  
  if ( grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
      .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)&
    call VecDestroy(iphase_all, ierr)

  if (grid%rk > 0.d0) then
    call VecDestroy(phis_all, ierr)
  endif
 ! if (iprint == 0) call VecScatterDestroy(scat_1dof, ierr)
 ! call VecScatterDestroy(scat_nph, ierr)

  close (IUNIT3)

  if (grid%iprint >= 1) then

!   if (ibrkcrv == 0) then

!   velocity fields
    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
                                vl_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
                             vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must 
! be followed by a call to VecDestroy(x_all)

#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(vl_nat, vl_all, ierr)
    call VecScatterCreateToZero(vl_nat, scat_3npdof, vl_all, ierr)
    call VecScatterBegin(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
#else
! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr)
    call VecScatterCreateToAll(vl_nat, scat_3npdof, vl_all, ierr)
    call VecScatterBegin(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                         SCATTER_FORWARD, ierr)
    call VecScatterEnd(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
    call VecScatterDestroy(scat_3npdof, ierr)
#endif

   ! endif
    
    if (grid%myrank == 0 .and. kplt > 0) then
       
      call VecGetArrayF90(vl_all, vl_p, ierr)
      
      ! write x,y,z-velocities
      ! note ordering of velocities: v(p,d,n) = v(p+(d-1)*Np+(n-1)*3*Np)
      if (kplt < 10) then
        write(fname,'(a10,i1,a4)') 'pflow_vel0', kplt, '.dat'
      else
        write(fname,'(a9,i2,a4)') 'pflow_vel', kplt, '.dat'
      endif
      write(*,*) '--> write output file: ',fname
      open(unit=IUNIT3,file=fname,action="write")
      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlx',q,'vgx',q,'vly',q,'vgy',q,'vlz',q, &
          'vgz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
            
            !x-connections
            nc = n
            nxm1 = nc-1
            ip = 0
            if (i>1 .and. i<grid%nx) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
            else if (i==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (i==grid%nx) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
            endif
            
            !y-connections
            nc = n
            nym1 = nc-grid%nx
            ip = 1
            if (j>1 .and. j<grid%ny) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            else if (j==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (j==grid%ny) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            endif
            
            !z-connections
            nc = n
            nzm1 = nc-grid%nxy
            ip = 2
            if (k>1 .and. k<grid%nz) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            else if (k==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (k==grid%nz) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            endif
            
            if(grid%nphase>1)then
              write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
              (vavel(ip),vaveg(ip),ip=1,3)
            else
              write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
              (vavel(ip),ip=1,3)
            endif 
!           write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
!           (vl_p(1+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
!           vl_p(2+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,ip=0,2)
          enddo
        enddo
      enddo
      close(IUNIT3)

!#if 0
      ! write x-velocities
      if (grid%nx > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlx0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlx', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' [''a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'y',q,'z',q,'vlx','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
   &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx-1,grid%ny,grid%nz
        do k = 1, grid%nz
          do j = 1, grid%ny
            do i = 1, grid%nx-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
               if(grid%nphase>1)then
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+(n-1)*3*grid%nphase)*grid%tconv,&
                 vl_p(2+(n-1)*3*grid%nphase)*grid%tconv
               else
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+(n-1)*3*grid%nphase)*grid%tconv!,&
                 !vl_p(2+(n-1)*3*grid%nphase)*grid%tconv
               endif                 
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif

      ! write y-velocities
      if (grid%ny > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vly0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vly', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vly','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny-1,grid%nz
        do k = 1, grid%nz
          do i = 1, grid%nx
            do j = 1, grid%ny-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
               if(grid%nphase>1)then
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
                 vl_p(2+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
               else
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
               endif                  
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif
      
      ! write z-velocities
      if (grid%nz > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlz0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlz', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlz','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz-1
        do j = 1, grid%ny
          do i = 1, grid%nx
            do k = 1, grid%nz-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
              if(grid%nphase>1)then
                write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                vl_p(1+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,&
                vl_p(2+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
              else
               write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                vl_p(1+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
              endif 
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif
!#endif

     call VecRestoreArrayF90(vl_all, vl_p, ierr)
    endif
   ! call VecRestoreArrayF90(vl_all, vl_p, ierr)
    call VecDestroy(vl_all, ierr)
  endif
  
! if (ibrkcrv >= 0 .or. iprint >= 1) then
!   call VecDestroy(vl_allMPI_Recv, ierr)
 !   call VecScatterDestroy(scat_3dof, ierr)
!   call VecDestroy(vl_nat, ierr)
! endif

  if (grid%iprint >= 2) call porperm_out (grid%t, grid%dt, grid%tconv, &
          kplt, grid%nx, grid%ny, grid%nz, grid%nmax, &
          grid%x, grid%y, grid%z, grid%flowsteps, scat_1dof, &
          grid%da_1_dof, grid%porosity, por_nat, grid%perm_xx, perm_nat, &
          grid%myrank)

  
    ! if (grid%iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)

  kplt = kplt + 1
  iplot = 0

  call VecDestroy(c_nat, ierr)
  call VecDestroy(phis_nat, ierr)
  call VecDestroy(t_nat, ierr)
  call VecDestroy(por_nat, ierr)
  call VecDestroy(perm_nat, ierr)
  call VecDestroy(p_nat, ierr)
  call VecDestroy(s_nat, ierr)
  call VecDestroy(vl_nat, ierr)
  call VecDestroy(x_nat, ierr)
  call VecDestroy(iphase_nat,ierr)
 
 ! if (iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)
!call VecScatterView(scat_1dof,PETSC_VIEWER_STDOUT_SELF)
end subroutine pflow_output

!======================================================================

  subroutine porperm_out(t, dt, tconv, kplt, nx, ny, nz, nmax, &
                         x, y, z, flowsteps, scat_1dof, da_1_dof, &
                         porosity, por_nat, perm, perm_nat, myrank)
  
  use PetscRelWrappers  ! For petsc-release compatibility.

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"
  
  integer, intent(in) :: kplt, nx, ny, nz, nmax, flowsteps, myrank
  real*8, intent(in) :: t, dt, tconv, x(:), y(:), z(:)
  Vec, intent(inout) :: porosity, por_nat, perm, perm_nat
  Vec :: por_all, perm_all
  
  DA :: da_1_dof
  
  VecScatter :: scat_1dof
  
  integer :: ierr, n
  
  real*8 :: fac, tyr
  
  real*8, pointer :: por_p(:), perm_p(:)
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab

#include "definitions.h"

  tab = char(9)
  q = '","'
  tyr =t/tconv

! write out porosity, permeability

! porosity field
  call DAGlobalToNaturalBegin(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)

! permeability field
! call DAGlobalToNaturalBegin(da_3, perm, INSERT_VALUES, perm_nat, ierr)
! call DAGlobalToNaturalEnd(da_3, perm, INSERT_VALUES, perm_nat, ierr)
  call DAGlobalToNaturalBegin(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)

! call VecScatterCreate(scatter, ierr)
#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(por_nat, por_all, ierr); CHKERRQ(ierr)
! call VecScatterCreateToZero(por_nat, scatter, por_all, ierr)
  call VecScatterBegin(scat_1dof, por_nat, por_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, por_nat, por_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
! call VecConvertMPIToMPIZero(perm_nat, perm_all, ierr); CHKERRQ(ierr)
  call VecScatterBegin(scat_1dof, perm_nat, perm_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, perm_nat, perm_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
#else
! call VecConvertMPIToSeqAll(por_nat, por_all, ierr); CHKERRQ(ierr)
! call VecScatterCreateToAll(por_nat, scatter, por_all, ierr)

! call VecScatterBegin(scat_1dof, por_nat, por_all, INSERT_VALUES, &
!                      SCATTER_FORWARD, ierr)
! call VecScatterEnd(scat_1dof, por_nat, por_all, INSERT_VALUES, &
!                    SCATTER_FORWARD, ierr)

! call VecConvertMPIToSeqAll(perm_nat, perm_all, ierr); CHKERRQ(ierr)

  call VecScatterCreateToAll(por_nat, scat_1dof, por_all, ierr)
  call VecScatterBegin(scat_1dof, por_nat, por_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, por_nat, por_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
! call VecScatterDestroy(scat_1dof, ierr)

  call VecScatterCreateToAll(perm_nat, scat_1dof, perm_all, ierr)
  call VecScatterBegin(scat_1dof, perm_nat, perm_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_1dof, perm_nat, perm_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
  call VecScatterDestroy(scat_1dof, ierr)
#endif

  if(myrank == 0) then
    
    if (kplt < 10) then
      write(fname,'(a10,i1,a4)') 'pflow_por0', kplt, '.dat'
    else
      write(fname,'(a9,i2,a4)') 'pflow_por', kplt, '.dat'
    endif
    
    write(*,*) '--> write output file: ',fname
    
    open(unit=IUNIT3,file=fname,action="write")
      
    call VecGetArrayF90(por_all, por_p, ierr)
    call VecGetArrayF90(perm_all, perm_p, ierr)
      
    fac = 1.d12
      
    if (nx+ny+nz .gt. nx*ny*nz) then ! 1D
        
      write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," y", &
 &    " dt =",1pe12.4," sec")') &
      flowsteps,t,tyr,dt
      write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
 &    "  -por- ","  -perm- ")')
!&    "  -por-      -permx-      -permy-      -permz-&
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (nz == 1) then ! 2D x-y

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,ny
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (ny == 1) then ! 2D x-z

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else ! 3D

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,nx,ny,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
    endif
    call VecRestoreArrayF90(por_all, por_p, ierr)
    call VecRestoreArrayF90(perm_all, perm_p, ierr)
    close (IUNIT3)
  endif
  call VecDestroy(por_all, ierr)
  call VecDestroy(perm_all, ierr)

  end subroutine porperm_out




 subroutine pflow_var_output(grid,kplt,iplot)
  
  use pflow_gridtype_module
  use TTPHASE_module
  use PetscRelWrappers  ! For petsc-release compatibility.

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"
!#include "include/mpiuni/mpif.h"
#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(inout) :: iplot, kplt
   
  Vec :: vl_nat, vl_all 
  real*8, pointer :: p_p(:), t_p(:), xx_p(:), phis_p(:), vl_p(:), var_p(:)
  real*8, pointer :: x_p(:), iphase_p(:)
  real*8, pointer :: fldflx(:), fldvol(:)
  real*8 xxx(1:grid%ndof)
  real*8, allocatable:: vvar(:), vavel(:, :)
  integer :: nxm1,nym1,nzm1, ndex, na, nv
  integer :: ierr,mm, nvar_out, status(MPI_STATUS_SIZE)
  integer :: i,ip,j,jn,k,ibrk,n,nc,nn, ipha, jj, ifound
  
  real*8 :: tyr, sum1, sum2, sum1v, sum2v, area, vol, vel, vf, iipha
  integer, save :: icall, icall_brk, icall_plt
   VecScatter :: scat_3npdof
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab
  
  character*6, allocatable :: var_name(:)!2+grid%nphase*(2+grid%nspec))
  
  data icall/1/, icall_brk/1/, icall_plt/1/
  
  nvar_out= 2+grid%nphase*(2+grid%nspec)
  allocate(var_name(nvar_out))
  
  var_name(1) = 'T'
  var_name(2) = 'p'
  do i = 1,grid%nphase
    write(var_name(i+2),'(a2,i3,a1)') 's(',i,')'
    write(var_name(grid%nphase + i+2), '(a2,i3,a1)') 'U(',i,')'
    do j =1, grid%nspec
      write(var_name(2*grid%nphase + (i-1)*grid%nphase + 2+j),'(a2,i1,a1,i1,a1)') 'X(',i,',',j,')'
    enddo
  enddo
  tab = char(9)
  q = '","'
  tyr = grid%t/grid%tconv

  size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
  size_var_node = size_var_use * (1+ grid%ndof)
  allocate(vvar(1:size_var_use))
  
  !save icall, icall_brk, icall_plt
  
! ibrkcrv = 0, no time-history
!         = 1, time-history output
  
! iplot  = 0, only time-history
!        = 1, spatial plot

! iprint = -1, none
!        =  0, p, t, s, c
!        =  1, vl
!        =  2, porosity, (permeability)
  !return
  
  if ((grid%ibrkcrv == 0 .and. iplot == 0) .or. grid%iprint == -1) then
    if (grid%iprint==-1 .and. iplot==1) then
      kplt = kplt + 1
      iplot = 0
    endif
    return
 endif
  
     
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1 .and. icall == 1) then
    icall = 0
    if (grid%myrank == 0) then
      write(fname,'(a9,a4)') 'pflow_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' iplot= ',iplot
      open(unit=IUNIT4,file=fname,action="write")
      write(IUNIT4,'("%#t            dt",100i12)') (i,i=1,grid%ibrkcrv), &
      (i,i=1,grid%ibrkcrv)
    endif
    if (iplot == 0) return
 endif
  
  if (grid%ndof == 1 .and. iplot == 0) return ! no time-history plot for 1 dof

! Create Natural Vec for output: use VecDuplicate here?
  
  ! time-history plots



! plot spatial data (iplot > 0)
  
    call VecGetArrayF90(grid%xx,  xx_p, ierr)
    call VecGetArrayF90(grid%var, var_p, ierr)
    call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  
  if (kplt < 10) then
      write(fname,'(a7,i1,a4)') 'pflow00', kplt, '.dat'
    else if (kplt < 100) then
      write(fname,'(a6,i2,a4)') 'pflow0', kplt, '.dat'
    else
      write(fname,'(a5,i3,a4)') 'pflow', kplt, '.dat'
    endif
    if(grid%myrank==0) write(*,*) '--> write output file: ',fname
    
    if(grid%myrank==0) open(unit=IUNIT3,file=fname,action="write")
      
    if (grid%nx+grid%ny+grid%nz .gt. grid%nx*grid%ny*grid%nz) then ! 1D
           if(grid%myrank ==0)then
        write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," [",a1,"]", &
     &    " dt =",1pe12.4," sec")') &
        grid%flowsteps,grid%t,tyr,grid%tunit,grid%dt
        if( grid%use_2ph == PETSC_TRUE)then
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
     &      "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
      &      "-xlco2-     -xgco2-        -Vf-     ")')
        elseif( grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
                .or. grid%use_flash == PETSC_TRUE &
                .or. grid%use_richards == PETSC_TRUE) then 
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   iphase ", &
     &      "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
      &      "-xlco2-     -xgco2-        -Vf-     ")')
        else
         write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        var_name(1), (q,var_name(i),i=2,nvar_out),'"'

        endif
      endif
    ! close(IUNIT3)

    ! call MPI_Bcast(IUNIT, 1, MPI_INTEGER, 0,PETSC_COMM_WORLD,ierr)
    ndex=1
    do na=0, grid%nmax-1
      ifound=0
      if(grid%myrank==0) then
        ifound=0
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            jn= 1 + (n-1)*grid%ndof 
            nv = 1 + (n-1) * size_var_node
            iipha= iphase_p(n)
            xxx= xx_p(jn:jn+grid%ndof-1)
            vvar= var_p(nv:nv + size_var_use-1)
            ndex=n
            exit
          endif
        enddo
      if(ifound==0)then  
        call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553,PETSC_COMM_WORLD, status,ierr)  
        call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
        na+na,PETSC_COMM_WORLD,status ,ierr)
      endif

    else
      do n=ndex,grid%nlmax
        if(grid%nL2A(n) == na) then
          jn= 1 + (n-1)*grid%ndof 
          nv = 1 + (n-1) * size_var_node
          iipha=int(iphase_p(n))
          call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
          xxx= xx_p(jn:jn+grid%ndof-1)
          call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
          call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
          MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
          ndex=n
          exit
        endif
      enddo
    endif
                
    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
    
    if(grid%myrank ==0 )then
      do i= 1, grid%nphase
        if(vvar(2+i)<=1D-30)then
          vvar(2+4*grid%nphase+i)=0.D0
          do j=1, grid%nspec
            vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
          enddo   
        endif
      enddo    
      write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
      vvar(1:2+grid%nphase), & ! Saturations
      vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
      vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
    endif
  enddo

  else if (grid%nz == 1) then ! 2D x-y
     if(grid%myrank ==0)then  
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                       .or. grid%use_vadose == PETSC_TRUE &
                                       .or. grid%use_flash == PETSC_TRUE &
                                        .or. grid%use_richards == PETSC_TRUE ) then
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'y',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
         else
     write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        var_name(1), (q,var_name(i),i=2,nvar_out),'"'
         endif
         write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
     &    '' , J='',i4)') tyr,grid%nx,grid%ny
         endif

        ndex=1
        do na=0, grid%nmax-1
          ifound=0
          if(grid%myrank==0) then
            ifound=0
            do n=ndex,grid%nlmax
              if(grid%nL2A(n) == na) then
                ifound=1
                jn= 1 + (n-1)*grid%ndof 
                nv = 1 + (n-1) * size_var_node
                iipha= iphase_p(n)
                xxx= xx_p(jn:jn+grid%ndof-1)
                vvar= var_p(nv:nv + size_var_use-1)
                ndex=n
              exit
            endif
          enddo
          if(ifound==0)then
            call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553, &
            PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
            na, PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
            MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr) 
          endif

            else
              do n=ndex,grid%nlmax
                if(grid%nL2A(n) == na) then
                  jn= 1 + (n-1)*grid%ndof 
                  nv = 1 + (n-1) * size_var_node
                  iipha=int(iphase_p(n))
                  call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
                  xxx= xx_p(jn:jn+grid%ndof-1)
                  call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
                  call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
                  MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
                  ndex=n
                  exit
                endif
              enddo  
            endif            

          call MPI_Barrier(PETSC_COMM_WORLD, ierr)

          if(grid%myrank ==0 )then
            do i= 1, grid%nphase
              if(vvar(2+i)<=1D-30)then
                vvar(2+4*grid%nphase+i)=0.D0
                do j=1, grid%nspec
                  vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
                enddo   
              endif
            enddo    
            write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
            vvar(1:2+grid%nphase), & ! Saturations
            vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
            vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
          endif
        enddo
          
    else if (grid%ny == 1) then ! 2D x-z

    if (grid%itecplot == 0) then
       if(grid%myrank==0) then
       write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      if( grid%use_2ph == PETSC_TRUE ) then
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
       elseif(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
            .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
      else
       write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        var_name(1), (q,var_name(i),i=2,nvar_out),'"'
      endif
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
   &      '' , K='',i4)') tyr,grid%nx,grid%nz
    endif

    ndex=1
    do na=0, grid%nmax-1
      ifound=0
      if(grid%myrank==0) then
        ifound=0
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            jn= 1 + (n-1)*grid%ndof 
            nv = 1 + (n-1) * size_var_node
            iipha= iphase_p(n)
            xxx= xx_p(jn:jn+grid%ndof-1)
            vvar= var_p(nv:nv + size_var_use-1)
            ndex=n
            exit
          endif
        enddo
        if(ifound==0)then  
          call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553, &
          PETSC_COMM_WORLD, status,ierr)  
          call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
          na, PETSC_COMM_WORLD, status,ierr)
          call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
          MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr) 
        endif

          else
            do n=ndex,grid%nlmax
              if(grid%nL2A(n) == na) then
                jn= 1 + (n-1)*grid%ndof 
                nv = 1 + (n-1) * size_var_node
                iipha=int(iphase_p(n))
                call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
                xxx= xx_p(jn:jn+grid%ndof-1)
                call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
                call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
                MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
                ndex=n
              exit
            endif
          enddo  
        endif

          call MPI_Barrier(PETSC_COMM_WORLD, ierr)

          if(grid%myrank ==0 )then
            do i= 1, grid%nphase
              if(vvar(2+i)<=1D-30)then
                vvar(2+4*grid%nphase+i)=0.D0
                do j=1, grid%nspec
                  vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
                enddo
              endif
            enddo
            write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
            vvar(1:2+grid%nphase), & ! Saturations
            vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
            vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
          endif
        enddo
      else
       
        if(grid%myrank==0) then
          write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
          if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   .or. grid%use_richards == PETSC_TRUE ) then
           write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'z',q,'phase',q,'pl',q,'pg',q,'T',q,'sl',q,'sg',q,'xl',q,'xg',q,'vf','"'
        else
          write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q, "ip",q,&
          var_name(1), (q,var_name(i),i=2,nvar_out),'"'
        endif
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &      '' , K='',i4)') tyr,grid%nx,grid%nz
        endif

        ndex=1
        do na=0, grid%nmax-1
          ifound=0
          if(grid%myrank==0) then
            ifound=0
            do n=ndex,grid%nlmax
              if(grid%nL2A(n) == na) then
                ifound=1
                jn= 1 + (n-1)*grid%ndof 
                nv = 1 + (n-1) * size_var_node
                iipha= iphase_p(n)
                xxx= xx_p(jn:jn+grid%ndof-1)
                vvar= var_p(nv:nv + size_var_use-1)
                ndex=n
                exit   
              endif
            enddo
            if(ifound==0)then  
              call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553, &
              PETSC_COMM_WORLD, status,ierr)
              call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
              na, PETSC_COMM_WORLD, status,ierr)
              call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
              MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr)
            endif
          else
            do n=ndex,grid%nlmax
              if(grid%nL2A(n) == na) then
                jn= 1 + (n-1)*grid%ndof 
                nv = 1 + (n-1) * size_var_node
                iipha=int(iphase_p(n))
                call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
                xxx= xx_p(jn:jn+grid%ndof-1)
                call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
                call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
                MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
                ndex=n
                exit
              endif
            enddo  
          endif

          call MPI_Barrier(PETSC_COMM_WORLD, ierr)

          if(grid%myrank ==0 )then
            do i= 1, grid%nphase
              if(vvar(2+i)<=1D-30)then
                vvar(2+4*grid%nphase+i)=0.D0
                do j=1, grid%nspec
                  vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
                enddo   
              endif
            enddo    
            write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
            vvar(1:2+grid%nphase), & ! Saturations
            vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
            vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
          endif
        enddo
      endif  ! endif of itechplot<>1

    else ! 3D
    if(grid%myrank==0)then
    write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
    if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   .or. grid%use_richards == PETSC_TRUE) then
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
    else
   write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        var_name(1), (q,var_name(i),i=2,nvar_out),'"'
     endif
    write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
    endif 


    ndex=1
    do na=0, grid%nmax-1
      ifound=0
      if(grid%myrank==0) then
        ifound=0
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            jn= 1 + (n-1)*grid%ndof 
            nv = 1 + (n-1) * size_var_node
            iipha= iphase_p(n)
            xxx= xx_p(jn:jn+grid%ndof-1)
            vvar= var_p(nv:nv + size_var_use-1)
            ndex=n
            exit
          endif
        enddo
      if(ifound==0)then  
        call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553,PETSC_COMM_WORLD, status,ierr)  
        call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr)
      endif

        else
          do n=ndex,grid%nlmax
            if(grid%nL2A(n) == na) then
              jn= 1 + (n-1)*grid%ndof 
              nv = 1 + (n-1) * size_var_node
              iipha=int(iphase_p(n))
              call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
              xxx= xx_p(jn:jn+grid%ndof-1)
              call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0, &
              na, PETSC_COMM_WORLD, ierr)  
              call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
              MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
              ndex=n
              exit
            endif
          enddo  
        endif            

        call MPI_Barrier(PETSC_COMM_WORLD, ierr)

        if(grid%myrank ==0 )then
          do i= 1, grid%nphase
            if(vvar(2+i)<=1D-30)then
              vvar(2+4*grid%nphase+i)=0.D0
              do j=1, grid%nspec
                vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
              enddo   
            endif
          enddo    
          write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
          vvar(1:2+grid%nphase), & ! Saturations
          vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
          vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
        endif
      enddo
    endif ! endif of geomtry selection

!   write out initial conditions

    if (grid%write_init == 1) then

      if(grid%myrank ==0)then
        fname = 'pflow_init0.dat'
        write(*,*) '--> write output file: ',fname
        close (IUNIT3)
        open(unit=IUNIT3,file=fname,action="write")
      if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   .or. grid%use_richards == PETSC_TRUE) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2", &
 &      "       p      ","      T      ","     sl(g)      ","       xl       xg   ")')
      else       
        write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        var_name(1), (q,var_name(i),i=2,nvar_out),'"'
      endif
    endif

    ndex=1
    do na=0, grid%nmax-1
      ifound=0
      if(grid%myrank==0) then
        ifound=0
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            jn= 1 + (n-1)*grid%ndof 
            nv = 1 + (n-1) * size_var_node
            iipha= iphase_p(n)
            xxx= xx_p(jn:jn+grid%ndof-1)
            vvar= var_p(nv:nv + size_var_use-1)
            ndex=n
            exit
          endif
        enddo
        if(ifound==0)then  
        call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553,PETSC_COMM_WORLD, status,ierr)  
        call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr) 
      endif

    else
      do n=ndex,grid%nlmax
        if(grid%nL2A(n) == na) then
          jn= 1 + (n-1)*grid%ndof 
          nv = 1 + (n-1) * size_var_node
          iipha=int(iphase_p(n))
          call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
          xxx= xx_p(jn:jn+grid%ndof-1)
          call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
          call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
          MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
          ndex=n
          exit
        endif
      enddo  
    endif

          call MPI_Barrier(PETSC_COMM_WORLD, ierr)

          if(grid%myrank ==0 )then
            do i= 1, grid%nphase
              if(vvar(2+i)<=1D-30)then
                vvar(2+4*grid%nphase+i)=0.D0
                do j=1, grid%nspec
                  vvar( 2+ 7*grid%nphase + (i-1)* grid%nspec + j)=0.D0
                enddo   
              endif
            enddo    
            write(IUNIT3,'(1p100e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), real(iipha), &
            vvar(1:2+grid%nphase), & ! Saturations
            vvar(3+4*grid%nphase:2+5*grid%nphase), &! Internal Energy
            vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec) !Mol fractions
          endif
        enddo
      endif
    
      
    
 
  ! if (iprint == 0) call VecScatterDestroy(scat_1dof, ierr)
 ! call VecScatterDestroy(scat_nph, ierr)

  close (IUNIT3)

    call VecRestoreArrayF90(grid%xx,  xx_p, ierr)
    call VecRestoreArrayF90(grid%var, var_p, ierr)
    call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)

  if (grid%iprint >= 1) then

!   if (ibrkcrv == 0) then

!   velocity fields

    call DACreateNaturalVector(grid%da_3np_dof, vl_nat, ierr)

    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
    vl_nat,ierr)
    call DAGlobalToNaturalEnd  (grid%da_3np_dof,grid%vl,INSERT_VALUES, &
    vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must 
! be followed by a call to VecDestroy(x_all)


#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(vl_nat, vl_all, ierr)
  call VecScatterCreateToZero(vl_nat, scat_3npdof, vl_all, ierr)
  call VecScatterBegin(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
#else
! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr)
  call VecScatterCreateToAll(vl_nat, scat_3npdof, vl_all, ierr)
  call VecScatterBegin(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                       SCATTER_FORWARD, ierr)
  call VecScatterEnd(scat_3npdof, vl_nat, vl_all, INSERT_VALUES, &
                     SCATTER_FORWARD, ierr)
 call VecScatterDestroy(scat_3npdof, ierr)
#endif




  ! endif
    
    if(grid%myrank == 0 .and. kplt > 0) then
       
      call VecGetArrayF90(vl_all, vl_p, ierr)
      allocate(vavel(3, grid%nphase))
      vavel=0.D0       
      ! write x,y,z-velocities
      ! note ordering of velocities: v(p,d,n) = v(p+(d-1)*Np+(n-1)*3*Np)
      if (kplt < 10) then
        write(fname,'(a10,i1,a4)') 'pflow_vel0', kplt, '.dat'
      else
        write(fname,'(a9,i2,a4)') 'pflow_vel', kplt, '.dat'
      endif
      write(*,*) '--> write output file: ',fname
      open(unit=IUNIT3,file=fname,action="write")
      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlx',q,'vgx',q,'vly',q,'vgy',q,'vlz',q, &
          'vgz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
                       !x-connections
            nc = n
            nxm1 = nc-1
            ip = 0
            do ipha=1, grid%nphase
              if (i>1 .and. i<grid%nx) then
                vavel(ip+1, ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
                + vl_p(ipha+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              else if (i==1) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              else if (i==grid%nx) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              endif
            enddo 

            !y-connections
            nc = n
            nym1 = nc-grid%nx
            ip = 1
            do ipha=1, grid%nphase
              if (j>1 .and. j<grid%ny) then
                vavel(ip+1,ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
                + vl_p(ipha+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              else if (j==1) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              else if (j==grid%ny) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              endif
            enddo

            !z-connections
            nc = n
            nzm1 = nc-grid%nxy
            ip = 2
            do ipha=1, grid%nphase
              if (k>1 .and. k<grid%nz) then
                vavel(ip+1,ipha) = 0.5d0*(vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase) &
                + vl_p(ipha+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              else if (k==1) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              else if (k==grid%nz) then
                vavel(ip+1,ipha) = (vl_p(ipha+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              endif
            enddo

            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
            ((vavel(ip,jj),jj=1,grid%nphase),ip=1,3)

!           write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
!           (vl_p(1+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
!           vl_p(2+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,ip=0,2)
          enddo
        enddo
      enddo
      close(IUNIT3)
      call VecRestoreArrayF90(vl_all, vl_p, ierr)
    endif
   ! call VecRestoreArrayF90(vl_all, vl_p, ierr)
    call VecDestroy(vl_all, ierr)
    call VecDestroy(vl_nat, ierr)
  endif
  
! if (ibrkcrv >= 0 .or. iprint >= 1) then
!   call VecDestroy(vl_all, ierr)
 !   call VecScatterDestroy(scat_3dof, ierr)
!   call VecDestroy(vl_nat, ierr)
! endif
 
  ! if (grid%iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)

  kplt = kplt + 1
  iplot = 0

! call VecDestroy(grid%c_nat, ierr)
! call VecDestroy(grid%phis_nat, ierr)
! call VecDestroy(grid%t_nat, ierr)
! call VecDestroy(grid%por_nat, ierr)
! call VecDestroy(grid%p_nat, ierr)
! call VecDestroy(grid%s_nat, ierr)
! call VecDestroy(grid%vl_nat, ierr)
! call VecDestroy(grid%x_nat, ierr)

  
  
 ! if (iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)
!call VecScatterView(scat_1dof,PETSC_VIEWER_STDOUT_SELF)
  end subroutine pflow_var_output

subroutine geh_io(grid, kplt)

  use pflow_gridtype_module
  use TTPHASE_module
  
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer :: kplt

  integer, save :: id = 0
  integer :: i, ierr
  character(len=20) :: id_string, filename
  Vec :: vec_1_dof
  PetscScalar, pointer :: vec_ptr(:)
  PetscViewer :: viewer

  if (id == 0) then
    open(unit=86,file='info.dat')
    if (associated(grid%imat)) then
      write(86,'(a8)') '+x +y +z'
    else
      write(86,'(a8)') '+x +y -z'
    endif
    write(86,'(i5,x,i5,x,i5)') grid%nx, grid%ny, grid%nz
  endif
  write(86,'(es20.8)') grid%t/(3600.d0*24.d0*365.d0)
  
  write(id_string,'(i7)') id
  
  call DACreateGlobalVector(grid%da_1_dof,vec_1_dof,ierr)

  if (kplt == 0) then
    ! x coordinates
    filename = 'xcoord.dat'
    call VecGetArrayF90(vec_1_dof, vec_ptr, ierr)
    do i=1,grid%nlmax
      vec_ptr(i) = grid%x(i)
    enddo
    call VecRestoreArrayF90(vec_1_dof, vec_ptr, ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
    call VecView(vec_1_dof,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    ! y coordinates
    filename = 'ycoord.dat'
    call VecGetArrayF90(vec_1_dof, vec_ptr, ierr)
    do i=1,grid%nlmax
      vec_ptr(i) = grid%y(i)
    enddo
    call VecRestoreArrayF90(vec_1_dof, vec_ptr, ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
    call VecView(vec_1_dof,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    ! z coordinates
    filename = 'zcoord.dat'
    call VecGetArrayF90(vec_1_dof, vec_ptr, ierr)
    do i=1,grid%nlmax
      vec_ptr(i) = grid%z(i)
    enddo
    call VecRestoreArrayF90(vec_1_dof, vec_ptr, ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
    call VecView(vec_1_dof,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  filename = 'pl_' // trim(adjustl(id_string)) // '.dat'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
  call VecStrideGather(grid%pressure,0,vec_1_dof,INSERT_VALUES,ierr)
  call VecView(vec_1_dof,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

!  filename = 'pg_' // trim(adjustl(id_string)) // '.dat'
!  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
!  call VecStrideGather(grid%pressure,1,vec_1_dof,INSERT_VALUES,ierr)
!  call VecView(vec_1_dof,viewer,ierr)
!  call PetscViewerDestroy(viewer,ierr)

  filename = 'sl_' // trim(adjustl(id_string)) // '.dat'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
  call VecStrideGather(grid%sat,0,vec_1_dof,INSERT_VALUES,ierr)
  call VecView(vec_1_dof,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

!  filename = 'sg_' // trim(adjustl(id_string)) // '.dat'
!  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
!  call VecStrideGather(grid%sat,1,vec_1_dof,INSERT_VALUES,ierr)
!  call VecView(vec_1_dof,viewer,ierr)
!  call PetscViewerDestroy(viewer,ierr)

  filename = 'xl_' // trim(adjustl(id_string)) // '.dat'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
  call VecStrideGather(grid%xmol,0,vec_1_dof,INSERT_VALUES,ierr)
  call VecView(vec_1_dof,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

!  filename = 'xg_' // trim(adjustl(id_string)) // '.dat'
!  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
!  call VecStrideGather(grid%xmol,1,vec_1_dof,INSERT_VALUES,ierr)
!  call VecView(vec_1_dof,viewer,ierr)
!  call PetscViewerDestroy(viewer,ierr)

  filename = 't_' // trim(adjustl(id_string)) // '.dat'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
  call VecView(grid%temp,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  filename = 'iph_' // trim(adjustl(id_string)) // '.dat'
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
  call VecView(grid%iphas,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  id = id + 1

  call VecDestroy(vec_1_dof,ierr)

end subroutine geh_io


subroutine pflow_soutput(grid,kplt,iplot)
  
  use pflow_gridtype_module
  use TTPHASE_module
  use PetscRelWrappers  ! For petsc-release compatibility.
  
  implicit none

#include "include/finclude/petsc.h"
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(inout) :: iplot, kplt
   
   
  Vec :: p_nat, t_nat, c_nat, phis_nat, por_nat, vl_nat, s_nat, x_nat, perm_nat,&
         iphase_nat

  real*4, pointer :: p_p(:), t_p(:), c_p(:), phis_p(:), vl_p(:), s_p(:)
  real*4, pointer :: x_p(:), iphase_p(:)
  real*8, pointer :: fldflx(:), fldvol(:)
  
  real*8 :: vavel(3),vaveg(3)
  integer :: nxm1,nym1,nzm1
  
  integer :: ierr
  integer :: i,ip,j,jn,k,ibrk,n,nc,nn
  integer :: vecsize
  
  real*8 :: tyr, sum1, sum2, sum1v, sum2v, area, vol, vel, vf
  integer, save :: icall, icall_brk, icall_plt
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab
  
  data icall/1/, icall_brk/1/, icall_plt/1/
  
   
  !save icall, icall_brk, icall_plt
  
! ibrkcrv = 0, no time-history
!         = 1, time-history output
  
! iplot  = -1, Glenn's format
! iplot  = 0, only time-history
!        = 1, spatial plot

! iprint = -1, none
!        =  0, p, t, s, c
!        =  1, vl
!        =  2, porosity, (permeability)
  !return

  if (iplot == 1 .and. grid%iprint == -2) then
    call geh_io(grid,kplt)
!    kplt = kplt + 1
!    iplot = 0
!    return
  endif
  
  if ((grid%ibrkcrv == 0 .and. iplot == 0) .or. grid%iprint == -1) then
    if (grid%iprint==-1 .and. iplot==1) then
      kplt = kplt + 1
      iplot = 0
    endif
    return
  endif
  
!  if(grid%nphase>1) call pflow_2phase_massbal(grid)
      
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1 .and. icall == 1) then
    icall = 0
    if (grid%myrank == 0) then
      write(fname,'(a9,a4)') 'pflow_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' iplot= ',iplot
      open(unit=IUNIT4,file=fname,action="write")
      write(IUNIT4,'("%#t            dt",100i12)') (i,i=1,grid%ibrkcrv), &
      (i,i=1,grid%ibrkcrv)
    endif
    if (iplot == 0) return
  endif
  
  if (grid%ndof == 1 .and. iplot == 0) return ! no time-history plot for 1 dof

! Create Natural Vec for output: use VecDuplicate here?
  call DACreateNaturalVector(grid%da_1_dof, c_nat, ierr)
  call VecDuplicate(c_nat, phis_nat, ierr)
  call VecDuplicate(c_nat, t_nat,    ierr)
  call VecDuplicate(c_nat, por_nat,  ierr)
  call VecDuplicate(c_nat, perm_nat,  ierr)
  call VecDuplicate(c_nat, iphase_nat, ierr)
  call DACreateNaturalVector(grid%da_nphase_dof, p_nat, ierr)
  call VecDuplicate(p_nat, s_nat, ierr)
  call VecDuplicate(p_nat, x_nat, ierr)
  call DACreateNaturalVector(grid%da_3np_dof, vl_nat, ierr)
  
  ! time-history plots
  
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1) then

!   concentration
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat, &
                                ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat,ierr)

!   velocity fields
    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat, &
                                ierr)
    call DAGlobalToNaturalEnd(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must be
! followed by a call to VecDestroy(x_all)

    call VecGetSize(c_nat, vecsize ,ierr)
    if (grid%myrank == 0) allocate(c_p(vecsize))
    call svectompizero(c_nat, c_p, ierr)

    call VecGetSize(vl_nat, vecsize, ierr)
    if (grid%myrank == 0) allocate(vl_p(vecsize))
    call svectompizero(vl_nat, vl_p, ierr)
   
    if (grid%myrank == 0) then
      allocate(fldflx(grid%ibrkcrv))
      allocate(fldvol(grid%ibrkcrv))
      
      fldflx = 0.d0
      fldvol = 0.d0
      do ibrk = 1, grid%ibrkcrv
        sum1 = 0.d0
        sum2 = 0.d0
        sum1v = 0.d0
        sum2v = 0.d0
        do k = grid%k1brk(ibrk),grid%k2brk(ibrk)
          do j = grid%j1brk(ibrk),grid%j2brk(ibrk)
            do i = grid%i1brk(ibrk),grid%i2brk(ibrk)
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
              vol = grid%dx0(i) * grid%dy0(j) * grid%dz0(k)
              if (grid%ibrkface(ibrk) == 1) then
                nn = i-1+(j-1)*grid%nx+(k-1)*grid%nxy
                area = grid%dy0(j) * grid%dz0(k)
              else if (grid%ibrkface(ibrk) == 2) then
                nn = i+(j-2)*grid%nx+(k-1)*grid%nxy
                area = grid%dx0(i)*grid%dz0(k)
              else if (grid%ibrkface(ibrk) == 3) then
                nn = i+(j-1)*grid%nx+(k-2)*grid%nxy
                area = grid%dx0(i)*grid%dy0(j)
              endif
              
              vel = vl_p(grid%ibrkface(ibrk)+3*(nn-1))

              sum1 = sum1 + vel*area*c_p(n)
              sum2 = sum2 + vel*area
              sum1v = sum1v + vol*c_p(n)
              sum2v = sum2v + vol
              
!             print *,'output-brk: ',grid%ibrkcrv,ibrk,i,j,k,n,nn, &
!             grid%ibrkface(ibrk),area, &
!             vol,vel*grid%tconv,c_p(n),sum1,sum2,sum1v,sum2v
            enddo
          enddo
        enddo
        if (sum2 .ne. 0.d0) fldflx(ibrk) = sum1/sum2
        if (sum2v .ne. 0.d0) fldvol(ibrk) = sum1v/sum2v
      enddo

      write(IUNIT4,'(1p100e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv, &
                                   (fldflx(i),i=1,grid%ibrkcrv), &
                                   (fldvol(i),i=1,grid%ibrkcrv)
      
      deallocate(c_p)
      deallocate(vl_p)
      deallocate(fldflx)
      deallocate(fldvol)
    endif

    if (iplot == 0)then
      call VecDestroy( c_nat, ierr)
      call VecDestroy( iphase_nat, ierr)
      call VecDestroy( phis_nat, ierr)
      call VecDestroy( t_nat,    ierr)
      call VecDestroy( por_nat,  ierr)
      call VecDestroy( p_nat, ierr)
      call VecDestroy( s_nat, ierr)
      call VecDestroy( x_nat, ierr)
      call VecDestroy( vl_nat, ierr)
      return
    endif
  endif
  
! plot spatial data (iplot > 0)
  
! presssure field
  call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%pressure,INSERT_VALUES, &
                              p_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%pressure,INSERT_VALUES, &
                            p_nat,ierr)

! temperature field
  call DAGlobalToNaturalBegin(grid%da_1_dof,grid%temp,INSERT_VALUES, &
                              t_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_1_dof,grid%temp,INSERT_VALUES, &
                            t_nat,ierr)

! saturation
  call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%sat,INSERT_VALUES, &
                              s_nat,ierr)
  call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%sat,INSERT_VALUES, &
                            s_nat,ierr)
  
! primary variables
  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE) then
    call DAGlobalToNaturalBegin(grid%da_nphase_dof,grid%xmol,INSERT_VALUES, &
                                x_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_nphase_dof,grid%xmol,INSERT_VALUES, &
                              x_nat,ierr)
  endif  
  
  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%iphas,INSERT_VALUES, &
                                iphase_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%iphas,INSERT_VALUES, &
                              iphase_nat,ierr)
  endif  
  
! concentration: check if already done in time-history plot
  if (grid%ibrkcrv == 0) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%conc,INSERT_VALUES, &
                                c_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%conc,INSERT_VALUES, &
                              c_nat,ierr)
  endif
  
! call VecView(phis,PETSC_VIEWER_STDOUT_WORLD,ierr)

! solid volume fraction
  if (grid%rk > 0.d0) then
    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%phis,INSERT_VALUES, &
                                phis_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%phis,INSERT_VALUES, &
                              phis_nat,ierr)
  endif
  
  icall_plt = 0

! note that the following routines call VecCreate(x_all) and therefore must be
! followed by a call to VecDestroy(x_all)

! call VecConvertMPIToSeqAll(p_nat, p_all, ierr);    CHKERRQ(ierr)
! call VecConvertMPIToSeqAll(t_nat, t_all, ierr);    CHKERRQ(ierr)
! call VecConvertMPIToSeqAll(s_nat, s_all, ierr);    CHKERRQ(ierr)

  call VecGetSize(p_nat, vecsize, ierr)
  if (grid%myrank == 0) allocate(p_p(vecsize))
  call svectompizero(p_nat, p_p, ierr)
  call VecDestroy(p_nat, ierr)

! call VecScatterCreateToAll(s_nat, scatter, s_all, ierr)
  call VecGetSize(s_nat, vecsize, ierr)
  if (grid%myrank == 0) allocate(s_p(vecsize))
  call svectompizero(s_nat, s_p, ierr)
  call VecDestroy(s_nat, ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
       .or. grid%use_richards == PETSC_TRUE) then
    call VecGetSize(x_nat, vecsize, ierr)
    if (grid%myrank == 0) allocate(x_p(vecsize))
    call svectompizero(x_nat, x_p, ierr)
    call VecDestroy(x_nat, ierr)
  endif
  
  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
      .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE) then
    call VecGetSize(iphase_nat, vecsize, ierr)
    if (grid%myrank == 0) allocate(iphase_p(vecsize))
    call svectompizero(iphase_nat, iphase_p, ierr)
    call VecDestroy(iphase_nat, ierr)
  endif

  call VecGetSize(t_nat, vecsize, ierr)
  if (grid%myrank == 0) allocate(t_p(vecsize))
  call svectompizero(t_nat, t_p, ierr)
  call VecDestroy(t_nat, ierr)

  call VecGetSize(c_nat, vecsize, ierr)
  if (grid%myrank == 0) allocate(c_p(vecsize))
  call svectompizero(c_nat, c_p, ierr)
  call VecDestroy(c_nat, ierr)

  if (grid%rk > 0.d0) then
    call VecGetSize(phis_nat, vecsize, ierr)
    if (grid%myrank == 0) allocate(phis_p(vecsize))
    call svectompizero(phis_nat, phis_p, ierr)
    call VecDestroy(phis_nat, ierr)
  endif

  tab = char(9)
  q = '","'
  tyr = grid%t/grid%tconv

  if (grid%myrank == 0) then
    if (grid%flowsteps == 0) then
      call SNESGetTolerances(grid%snes, grid%atol, grid%rtol, grid%stol, &
                             grid%maxit, grid%maxf, ierr)
!       write(*,'("%#atol, rtol, stol= ",1p3e12.4," maxit, maxf= ",2i6)') &
!       grid%atol,grid%rtol,grid%stol,grid%maxit,grid%maxf
  
!     Calculate the x, y, z vectors that give the 
!     physical coordinates of each cell.
!     call pflowGrid_compute_xyz(grid)
    endif
    
    if (kplt < 10) then
      write(fname,'(a7,i1,a4)') 'pflow00', kplt, '.dat'
    else if (kplt < 100) then
      write(fname,'(a6,i2,a4)') 'pflow0', kplt, '.dat'
    else
      write(fname,'(a5,i3,a4)') 'pflow', kplt, '.dat'
    endif
    write(*,*) '--> write output file: ',fname
    open(unit=IUNIT3,file=fname,action="write")
      
    if (grid%nx+grid%ny+grid%nz .gt. grid%nx*grid%ny*grid%nz) then ! 1D
    
      write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," [",a1, &
        & "]", " dt =",1pe12.4," sec")') &
        grid%flowsteps,grid%t,tyr,grid%tunit,grid%dt
      if (grid%use_2ph == PETSC_TRUE) then
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
          & "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
          & "-xlco2-     -xgco2-        -Vf-     ")')
      elseif( grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
              .or. grid%use_flash == PETSC_TRUE) then 
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   iphase ", &
          & "  -pl[Pa]-     -pg[Pa]-     -t[C]-      -sl-        -sg-    ", &
          & "-xlco2-     -xgco2-        -Vf-     ")')
      else
        write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
          & "  -p[Pa]-      -t[C]-      -sl(g)-      -C[mol/L]-      -Vf- &
          & ")')
      endif

      do n = 1, grid%nmax
        jn = grid%jh2o+(n-1)*grid%nphase
        vf = 0.d0
        if (grid%rk > 0.d0) vf = phis_p(n)
        if (grid%use_2ph == PETSC_TRUE ) then
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn),s_p(jn+1), &
            x_p(jn),x_p(jn+1),vf
        elseif (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
                .or. grid%use_flash == PETSC_TRUE)then !.or. grid%use_richards == PETSC_TRUE)then
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            iphase_p(n),&
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn:jn+grid%nphase-1), &
            x_p(jn:jn+grid%nphase-1),vf
      !    print *,'output', grid%nphase, n, (p_p(j),j=jn,jn+grid%nphase-grid%jh2o),t_p(n),s_p(jn:jn+grid%nphase-1), &
      !      x_p(jn:jn+grid%nphase-1),vf
        else
          write(IUNIT3,'(1p100e12.4)') grid%x(n), grid%y(n), grid%z(n), &
            (p_p(j),j=jn,jn+grid%nphase-grid%jh2o), t_p(n), s_p(jn), c_p(n), vf
        endif
      enddo

    else if (grid%nz == 1) then ! 2D x-y

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
          grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE)then
  
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
      else
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
      endif
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
        & '' , J='',i4)') tyr,grid%nx,grid%ny

      do n = 1, grid%nmax
        jn = grid%jh2o+(n-1)*grid%nphase
        vf = 0.d0
        if (grid%rk > 0.d0) vf = phis_p(n)
        if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
            grid%use_vadose == PETSC_TRUE.or.grid%use_flash == PETSC_TRUE &
             ) then
          write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n), &
             p_p(jn+1), t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
        else
          write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), &
            p_p(jn), t_p(n), s_p(jn), c_p(n), vf
        endif
      enddo
          
    else if (grid%ny == 1) then ! 2D x-z

      if (grid%itecplot == 0) then
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        if( grid%use_2ph == PETSC_TRUE ) then
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
         elseif(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
                .or. grid%use_flash == PETSC_TRUE )then
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
        else
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
        endif
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
          & '' , K='',i4)') tyr,grid%nx,grid%nz
        
        do n = 1, grid%nmax
          jn = grid%jh2o+(n-1)*grid%nphase
          vf = 0.d0
          if (grid%rk > 0.d0) vf = phis_p(n)
          if ( grid%use_2ph == PETSC_TRUE) then
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn+1), t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
          elseif(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
                 .or. grid%use_flash == PETSC_TRUE &
                 ) then
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n),iphase_p(n), &
              p_p(jn+1),t_p(n), s_p(jn), x_p(jn), x_p(jn+1),vf
          else
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn), t_p(n), s_p(jn), c_p(n), vf
          endif
        enddo

      else

        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
            grid%use_vadose == PETSC_TRUE.or.grid%use_flash == PETSC_TRUE &
             ) then
           write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
             'x',q,'z',q,'phase',q,'pl',q,'pg',q,'T',q,'sl',q,'sg',q,'xl',q, &
             'xg',q,'vf','"'
        else
          write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
        endif
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
          & '' , J='',i4)') tyr,grid%nx,grid%nz
        
!        do i = 1, grid%nx
!          do k = 1, grid%nz
!            n = i+(k-1)*grid%nxy
        do n=1, grid%nmax
          jn = grid%jh2o+(n-1)*grid%nphase
          vf = 0.d0
          if (grid%rk > 0.d0) vf = phis_p(n)
            
!            print *,'pflow_output: ',n,jn,p_p(jn),p_p(jn+1), t_p(n), &
!                                    s_p(jn),s_p(jn+1), x_p(jn),x_p(jn+1)
            
          if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
              .or. grid%use_vadose == PETSC_TRUE &
              .or. grid%use_flash == PETSC_TRUE &
              ) then
            write(IUNIT3,'(1p100e12.4)') grid%x(n),grid%z(n),iphase_p(n), &
              p_p(jn),p_p(jn+1),t_p(n),s_p(jn),s_p(jn+1),x_p(jn),x_p(jn+1),vf
          else
            write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%z(n), &
              p_p(jn), t_p(n), s_p(jn), c_p(n), vf
          endif
!          enddo
        enddo
      endif

    else ! 3D

    write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
    if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   ) then
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'y',q,'z',q,'phase',q,'p',q,'T',q,'sl(g)',q,'xl',q,'xg',q,'vf','"'
    else
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
        'x',q,'y',q,'z',q,'p',q,'T',q,'sl(g)',q,'c',q,'vf','"'
    endif
    write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
      & '' , J='',i4,'' K='',i4)') tyr,grid%nx,grid%ny,grid%nz
    do n = 1,  grid%nmax
      jn = grid%jh2o+(n-1)*grid%nphase
      vf = 0.d0
      if (grid%rk > 0.d0) vf = phis_p(n)
      if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
          .or. grid%use_flash == PETSC_TRUE &   
          .or. grid%use_vadose == PETSC_TRUE) then
        write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), grid%z(n), &
          iphase_p(n), p_p(jn+1), t_p(n), s_p(jn+1), x_p(jn),x_p(jn+1), vf
      else
        write(IUNIT3,'(1p10e12.4)') grid%x(n), grid%y(n), grid%z(n), &
          p_p(jn), t_p(n), s_p(jn), c_p(n), vf
      endif
    enddo
  endif

!   write out initial conditions

  if (grid%write_init == 1) then
    fname = 'pflow_init0.dat'
    write(*,*) '--> write output file: ',fname
    close (IUNIT3)
    open(unit=IUNIT3,file=fname,action="write")
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
        .or. grid%use_flash == PETSC_TRUE & 
        .or. grid%use_vadose == PETSC_TRUE &
        .or. grid%use_richards == PETSC_TRUE) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2", &
          & "       p      ","      T      ","     sl(g)      ", &
          & "       xl       xg   ")')
      else       
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2", &
          & "       p      ","      T      ","     sl(g)      ", &
          & "      C      ")')
      endif
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
            jn = grid%jh2o+(n-1)*grid%nphase
            if (grid%use_2ph == PETSC_TRUE .or. &
                grid%use_mph == PETSC_TRUE .or. &
                grid%use_flash == PETSC_TRUE .or. &
                grid%use_vadose == PETSC_TRUE .or.&
                grid%use_richards == PETSC_TRUE ) then
              write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') i,i,j,j,k,k, &
                p_p(jn), t_p(n), s_p(jn), x_p(jn),x_p(jn+1)
            else
              write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') i,i,j,j,k,k, &
                p_p(jn), t_p(n), s_p(jn), c_p(n)
            endif
          enddo
        enddo
      enddo
      write(IUNIT3,'("/",/)')
    endif
    
    deallocate(p_p)
    deallocate(c_p)
    deallocate(s_p)
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
        .or. grid%use_flash == PETSC_TRUE &
        .or. grid%use_vadose == PETSC_TRUE &
         .or. grid%use_richards == PETSC_TRUE ) then
      deallocate(x_p)
      deallocate(iphase_p)
    endif
    if (grid%rk > 0.d0) then
      deallocate(phis_p)
    endif
  endif
  
  close (IUNIT3)

  if (grid%iprint >= 1) then

!   if (ibrkcrv == 0) then

!   velocity fields
    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
                                vl_nat,ierr)
    call DAGlobalToNaturalEnd(grid%da_3np_dof,grid%vl,INSERT_VALUES, &
                             vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must 
! be followed by a call to VecDestroy(x_all)

! call VecConvertMPIToSeqAll(vl_nat, vl_all, ierr)
    call VecGetSize(vl_nat, vecsize, ierr)
    if (grid%myrank == 0) allocate(vl_p(vecsize))
    call svectompizero(vl_nat, vl_p, ierr)
    call VecDestroy(vl_nat, ierr)

   ! endif
    
    if (grid%myrank == 0 .and. kplt > 0) then
       
      ! write x,y,z-velocities
      ! note ordering of velocities: v(p,d,n) = v(p+(d-1)*Np+(n-1)*3*Np)
      if (kplt < 10) then
        write(fname,'(a10,i1,a4)') 'pflow_vel0', kplt, '.dat'
      else
        write(fname,'(a9,i2,a4)') 'pflow_vel', kplt, '.dat'
      endif
      write(*,*) '--> write output file: ',fname
      open(unit=IUNIT3,file=fname,action="write")
      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlx',q,'vgx',q,'vly',q,'vgy',q,'vlz',q, &
          'vgz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            n = i+(j-1)*grid%nx+(k-1)*grid%nxy
            
            !x-connections
            nc = n
            nxm1 = nc-1
            ip = 0
            if (i>1 .and. i<grid%nx) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
            else if (i==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (i==grid%nx) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nxm1-1)*3*grid%nphase))*grid%tconv
            endif
            
            !y-connections
            nc = n
            nym1 = nc-grid%nx
            ip = 1
            if (j>1 .and. j<grid%ny) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            else if (j==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (j==grid%ny) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nym1-1)*3*grid%nphase))*grid%tconv
            endif
            
            !z-connections
            nc = n
            nzm1 = nc-grid%nxy
            ip = 2
            if (k>1 .and. k<grid%nz) then
              vavel(ip+1) = 0.5d0*(vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(1+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = 0.5d0*(vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase) &
              + vl_p(2+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            else if (k==1) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nc-1)*3*grid%nphase))*grid%tconv
            else if (k==grid%nz) then
              vavel(ip+1) = (vl_p(1+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
              if(grid%nphase>1) vaveg(ip+1) = (vl_p(2+ip*grid%nphase+(nzm1-1)*3*grid%nphase))*grid%tconv
            endif
            
            if(grid%nphase>1)then
              write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
              (vavel(ip),vaveg(ip),ip=1,3)
            else
              write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
              (vavel(ip),ip=1,3)
            endif 
!           write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
!           (vl_p(1+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
!           vl_p(2+ip*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,ip=0,2)
          enddo
        enddo
      enddo
      close(IUNIT3)

      ! write x-velocities
      if (grid%nx > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlx0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlx', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' [''a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'y',q,'z',q,'vlx','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
   &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx-1,grid%ny,grid%nz
        do k = 1, grid%nz
          do j = 1, grid%ny
            do i = 1, grid%nx-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
               if(grid%nphase>1)then
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+(n-1)*3*grid%nphase)*grid%tconv,&
                 vl_p(2+(n-1)*3*grid%nphase)*grid%tconv
               else
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+(n-1)*3*grid%nphase)*grid%tconv!,&
                 !vl_p(2+(n-1)*3*grid%nphase)*grid%tconv
               endif                 
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif

      ! write y-velocities
      if (grid%ny > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vly0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vly', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vly','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny-1,grid%nz
        do k = 1, grid%nz
          do i = 1, grid%nx
            do j = 1, grid%ny-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
               if(grid%nphase>1)then
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv, &
                 vl_p(2+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
               else
                 write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                 vl_p(1+grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
               endif                  
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif
      
      ! write z-velocities
      if (grid%nz > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlz0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlz', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=IUNIT3,file=fname,action="write")
        write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlz','"'
        write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz-1
        do j = 1, grid%ny
          do i = 1, grid%nx
            do k = 1, grid%nz-1
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy
              if(grid%nphase>1)then
                write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                vl_p(1+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv,&
                vl_p(2+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
              else
               write(IUNIT3,'(1p10e12.4)') grid%x(n),grid%y(n),grid%z(n), &
                vl_p(1+2*grid%nphase+(n-1)*3*grid%nphase)*grid%tconv
              endif 
            enddo
          enddo
        enddo
        close(IUNIT3)
      endif

     deallocate(vl_p)
    endif
  endif
  
! if (ibrkcrv >= 0 .or. iprint >= 1) then
!   call VecDestroy(vl_allMPI_Recv, ierr)
 !   call VecScatterDestroy(scat_3dof, ierr)
!   call VecDestroy(vl_nat, ierr)
! endif

  if (grid%iprint >= 2) call porperm_sout (grid%t, grid%dt, grid%tconv, &
          kplt, grid%nx, grid%ny, grid%nz, grid%nmax, &
          grid%x, grid%y, grid%z, grid%flowsteps, &
          grid%da_1_dof, grid%porosity, por_nat, grid%perm_xx, perm_nat, &
          grid%myrank)

  
    ! if (grid%iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)

  kplt = kplt + 1
  iplot = 0

  call VecDestroy(por_nat, ierr)
  call VecDestroy(perm_nat, ierr)
 
 ! if (iprint >= 1) call VecScatterDestroy(scat_1dof, ierr)
!call VecScatterView(scat_1dof,PETSC_VIEWER_STDOUT_SELF)
end subroutine pflow_soutput


!======================================================================
! This is the single-precision version of porperm_out().  --RTM
  subroutine porperm_sout(t, dt, tconv, kplt, nx, ny, nz, nmax, &
                         x, y, z, flowsteps, da_1_dof, &
                         porosity, por_nat, perm, perm_nat, myrank)
  
  use PetscRelWrappers  ! For petsc-release compatibility.

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

#include "definitions.h"
  
  integer, intent(in) :: kplt, nx, ny, nz, nmax, flowsteps, myrank
  real*8, intent(in) :: t, dt, tconv, x(:), y(:), z(:)
  Vec, intent(inout) :: porosity, por_nat, perm, perm_nat
  
  DA :: da_1_dof
  
  integer :: ierr, n
  integer :: vecsize
  
  real*8 :: fac, tyr
  
  real*4, pointer :: por_p(:), perm_p(:)
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab

#include "definitions.h"

  tab = char(9)
  q = '","'
  tyr =t/tconv

! write out porosity, permeability

! porosity field
  call DAGlobalToNaturalBegin(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, porosity, INSERT_VALUES, &
              por_nat, ierr)

! permeability field
! call DAGlobalToNaturalBegin(da_3, perm, INSERT_VALUES, perm_nat, ierr)
! call DAGlobalToNaturalEnd(da_3, perm, INSERT_VALUES, perm_nat, ierr)
  call DAGlobalToNaturalBegin(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)
  call DAGlobalToNaturalEnd(da_1_dof, perm, INSERT_VALUES, &
              perm_nat, ierr)

! call VecScatterCreate(scatter, ierr)
! call VecConvertMPIToSeqAll(por_nat, por_all, ierr); CHKERRQ(ierr)
! call VecScatterCreateToAll(por_nat, scatter, por_all, ierr)

! call VecScatterBegin(scat_1dof, por_nat, por_all, INSERT_VALUES, &
!                      SCATTER_FORWARD, ierr)
! call VecScatterEnd(scat_1dof, por_nat, por_all, INSERT_VALUES, &
!                    SCATTER_FORWARD, ierr)

! call VecConvertMPIToSeqAll(perm_nat, perm_all, ierr); CHKERRQ(ierr)

  call VecGetSize(por_nat, vecsize, ierr)
  if (myrank == 0) allocate(por_p(vecsize))
  call svectompizero(por_nat, por_p, ierr)
  call VecDestroy(por_nat, ierr)

  call VecGetSize(perm_nat, vecsize, ierr)
  if (myrank == 0) allocate(perm_p(vecsize))
  call svectompizero(perm_nat, perm_p, ierr)
  call VecDestroy(perm_nat, ierr)

  if(myrank == 0) then
    
    if (kplt < 10) then
      write(fname,'(a10,i1,a4)') 'pflow_por0', kplt, '.dat'
    else
      write(fname,'(a9,i2,a4)') 'pflow_por', kplt, '.dat'
    endif
    
    write(*,*) '--> write output file: ',fname
    
    open(unit=IUNIT3,file=fname,action="write")
      
    fac = 1.d12
      
    if (nx+ny+nz .gt. nx*ny*nz) then ! 1D
        
      write(IUNIT3,'("%#step: ",i6," time=",1pe12.4," sec,",1pe12.4," y", &
 &    " dt =",1pe12.4," sec")') &
      flowsteps,t,tyr,dt
      write(IUNIT3,'("%#   -x[m]-      -y[m]-      -z[m]-   ", &
 &    "  -por- ","  -perm- ")')
!&    "  -por-      -permx-      -permy-      -permz-&
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (nz == 1) then ! 2D x-y

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,ny
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else if (ny == 1) then ! 2D x-z

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4)') tyr,nx,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
          
    else ! 3D

      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,nx,ny,nz
        
      do n = 1, nmax
        write(IUNIT3,'(1p10e12.4)') x(n), y(n), z(n), &
        por_p(n),perm_p(n) !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
      enddo
    endif
    deallocate(por_p)
    deallocate(perm_p)
    close (IUNIT3)
  endif

  end subroutine porperm_sout
  
end module pflow_output_module
