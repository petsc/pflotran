
module pflow_output_module_new
  use pflow_chkptheader

private
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

#if (PETSC_VERSION_RELEASE == 0 || PETSC_VERSION_SUBMINOR == 3)
#include "include/finclude/petscbag.h"

  Interface PetscBagGetData
    Subroutine PetscBagGetData(bag,ctx,ierr)
      use pflow_chkptheader
      PetscBag bag
      type(pflowChkPtHeader), pointer :: ctx
      PetscErrorCode ierr
    End Subroutine
  End Interface PetscBagGetData
#endif


  Vec :: var_plot, vl_plot_loc
  real*8, pointer :: vvar(:),vl_plot_loc_p(:) 
  integer, private, save :: size_var_use, size_var_node
  
 public   pflow_output_new
  
  contains
  
  subroutine pflow_output_new(grid, kplt,iplot)
 
   use pflow_gridtype_module
  use TTPHASE_module
  use PetscRelWrappers  ! For petsc-release compatibility.
  use pflow_output_module, only: geh_io
  use pflow_checkpoint
  implicit none


#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  integer, intent(inout) :: iplot, kplt
  PetscViewer viewer
#if (PETSC_VERSION_RELEASE == 0 || PETSC_VERSION_SUBMINOR == 3)
  PetscBag bag
  type(pflowChkPtHeader), pointer :: header
#endif

  real*8, pointer :: p_p(:), t_p(:), c_p(:), phis_p(:), vl_p(:), s_p(:)
  real*8, pointer :: x_p(:), iphase_p(:), xx_p(:), var_p(:), var_plot_p(:)
  real*8, pointer :: fldflx(:), fldvol(:)
  real*8 pres(grid%nphase), temp, conc, xmol(grid%nphase), sat(grid%nphase)
  
  real*8 :: vavel(3*grid%nphase), velo(3*grid%nphase)
  integer :: nxm1,nym1,nzm1
  
  integer :: ierr, ifound, status(MPI_STATUS_SIZE)
  integer :: i,ip,j,jn,k,ibrk,n,nc,nn, nz,ndex, ndex2, na
  integer :: nv, iipha, ng, ng_upstream
  
  real*8 :: tyr, sum1, sum2, sum1v, sum2v, area, vol, vel, vf
  integer, save :: icall, icall_brk, icall_plt
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab

  real*8 xxx(1:grid%ndof)

  integer ntstep, iflgcut, ihalcnt, its
  
  PetscTruth :: use_soutput
  DA ::  da_var_plot
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

 ! call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_soutput", &
 !                          use_soutput, ierr)
 ! if (use_soutput == PETSC_TRUE) then
 !   call pflow_soutput(grid, kplt, iplot)
 !   return
 ! endif
  
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
      
  if(icall==1)then
    call DACreateLocalVector(grid%da_3np_dof, vl_plot_loc, ierr)
    if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
        .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then 
    
        size_var_use = (2+7*grid%nphase + 2* grid%nspec*grid%nphase)
        size_var_node = size_var_use *(grid%ndof +1)
          call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    grid%nx,grid%ny,grid%nz,grid%npx,grid%npy,grid%npz,&
                    size_var_use,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    da_var_plot,ierr)
          call DACreateGlobalVector(da_var_plot, var_plot, ierr)
         ! print *, ' created varplot'
           if (grid%itecplot ==PETSC_TRUE)then
             allocate(vvar(1:size_var_use))
            
       !      if(grid%iprint>=1)allocate(vavel(3*grid%nphase))
          endif
     endif         
  endif     



      
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1 .and. icall == 1) then
    
    if (grid%myrank == 0) then
      write(fname,'(a9,a4)') 'pflow_his','.dat'
      write(*,*) '--> open time-history file: ',fname,' iplot= ',iplot
      open(unit=IUNIT4,file=fname,action="write")
      write(IUNIT4,'("%#t            dt",100i12)') (i,i=1,grid%ibrkcrv), &
      (i,i=1,grid%ibrkcrv)
    endif
    if (iplot == 0) return
  endif
  
  icall = 0
  if (grid%ndof == 1 .and. iplot == 0) return ! no time-history plot for 1 dof

    
    
! time-history plots
  
  if (grid%ibrkcrv > 0 .and. grid%ndof > 1) then
  print *, 'begin his' 
!   concentration
!    call DAGlobalToNaturalBegin(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat, &
!                                ierr)
!    call DAGlobalToNaturalEnd(grid%da_1_dof,grid%conc,INSERT_VALUES,c_nat,ierr)

!   velocity fields
!    call DAGlobalToNaturalBegin(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat, &
 !                               ierr)
 !   call DAGlobalToNaturalEnd(grid%da_3np_dof,grid%vl,INSERT_VALUES,vl_nat,ierr)

! note that the following routines call VecCreate(x_all) and therefore must be
! followed by a call to VecDestroy(x_all)


      call VecGetArrayF90(grid%conc, c_p, ierr)
      call VecGetArrayF90(grid%vl, vl_p, ierr)
      
    ! print *,'pflow_output_new: ',vl_p
      
      if(grid%myrank == 0)then
        allocate(fldflx(grid%ibrkcrv))
        allocate(fldvol(grid%ibrkcrv))
        fldflx = 0.d0
        fldvol = 0.d0
      endif  
      
        do ibrk = 1, grid%ibrkcrv
        sum1 = 0.d0
        sum2 = 0.d0
        sum1v = 0.d0
        sum2v = 0.d0
        ndex=1; ndex2 = 1
        do k = grid%k1brk(ibrk),grid%k2brk(ibrk)
          do j = grid%j1brk(ibrk),grid%j2brk(ibrk)
            do i = grid%i1brk(ibrk),grid%i2brk(ibrk)
              na = i+(j-1)*grid%nx+(k-1)*grid%nxy-1
              
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
              if (nn==0) nn=1
! look up for velocity
              ! print *, 'pflow_output_new: vel', nn        
              if(grid%myrank==0) then
                  ifound=0 
                  
                  do n=ndex,grid%nlmax
                    if(grid%nL2A(n) == (nn-1)) then
                      ifound=1
                      vel = vl_p(1+ (grid%ibrkface(ibrk)-1)*grid%nphase + 3*grid%nphase*(n-1))
                      ndex=n
                      
            !         print *,'pflow_output_new: vel',n,nn,ifound,grid%ibrkface(ibrk),vel,&
            !           vl_p(1+3*(n-1):3*n)
                      exit
                    endif
                  enddo
                  if(ifound==0)then  
                    call MPI_Recv(vel,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                                  na, PETSC_COMM_WORLD, status,ierr)
                   endif
               else
                     do n=ndex,grid%nlmax
                       if(grid%nL2A(n) == (nn-1)) then
                          ifound=1
                          vel = vl_p(1+ (grid%ibrkface(ibrk)-1)*grid%nphase + 3*grid%nphase*(n-1))                        
                           
                          call MPI_Send(vel, 1, MPI_DOUBLE_PRECISION,0, &
                             na, PETSC_COMM_WORLD, ierr)  

                          ndex=n
                          exit
                       endif
                     enddo
               endif

!look up for concentration

              if(grid%myrank==0) then
                  ifound=0
                  do n=ndex2,grid%nlmax
                    if(grid%nL2A(n) == na) then
                      ifound=1
                      conc=  c_p(n)
                      ndex2=n
                      exit
                    endif
                  enddo
                  if(ifound==0)then  
                    call MPI_Recv(conc, 1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                                  na, PETSC_COMM_WORLD, status,ierr)
                    endif
               else
                     do n=ndex2,grid%nlmax
                       if(grid%nL2A(n) == na) then
                          ifound=1
                          conc=  c_p(n)
                          call MPI_Send(conc, 1, MPI_DOUBLE_PRECISION,0, &
                             na, PETSC_COMM_WORLD, ierr)  
                        
                          ndex2=n
                          exit
                       endif
                     enddo
               endif


                
    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
              
           if(grid%myrank ==0)then
              sum1 = sum1 + vel*area*conc
              sum2 = sum2 + vel*area
              sum1v = sum1v + vol*conc
              sum2v = sum2v + vol
           endif   
!             print *,'output-brk: ',grid%nphase,grid%ibrkcrv,ibrk,i,j,k,n,nn, &
!             grid%ibrkface(ibrk),area, &
!             vol,vel*grid%tconv,conc,sum1,sum2,sum1v,sum2v
            enddo
          enddo
        enddo
       if(grid%myrank ==0)then
        if (sum2 .ne. 0.d0) fldflx(ibrk) = sum1/sum2
        if (sum2v .ne. 0.d0) fldvol(ibrk) = sum1v/sum2v
       endif 
    enddo ! end loop of ibrkcrv

    if(grid%myrank ==0)then
      write(IUNIT4,'(1p100e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv, &
                                   (fldflx(i),i=1,grid%ibrkcrv), &
                                   (fldvol(i),i=1,grid%ibrkcrv) 
      
      deallocate(fldflx)
      deallocate(fldvol)
    endif
    
      call VecRestoreArrayF90(grid%conc, c_p, ierr)
      call VecRestoreArrayF90(grid%vl, vl_p, ierr)
      

    
    if (iplot == 0)then
        return
    endif
  endif


! plot  spatial data (iplot > 0)

 icall_plt = 0

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
  endif

! Varibale Get pointers

!print *, 'begin rebuild output var'
 if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
        .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
    call VecGetArrayF90(grid%var, var_p, ierr)
    call VecGetArrayF90(var_plot, var_plot_p, ierr)
    do n= 1, grid%nlmax
      var_plot_p(1+(n-1)*size_var_use: n*size_var_use) = var_p(1+(n-1)*size_var_node: (n-1)* size_var_node+size_var_use)
    enddo   
    call VecRestoreArrayF90(grid%var, var_p, ierr)
    call VecRestoreArrayF90(var_plot, var_plot_p, ierr)
  endif

! print *, 'end rebuild output var'

 if (grid%itecplot ==PETSC_TRUE)then
   if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
        .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
      call VecGetArrayF90(grid%iphas, iphase_p, ierr)
      call VecGetArrayF90(grid%xx, xx_p, ierr)
      call VecGetArrayF90(grid%var, var_p, ierr)
    else
      call VecGetArrayF90(grid%pressure, p_p, ierr)
      call VecGetArrayF90(grid%temp, t_p, ierr)
      call VecGetArrayF90(grid%conc, c_p, ierr)
      call VecGetArrayF90(grid%sat, s_p, ierr)
      if (grid%use_2ph == PETSC_TRUE) call VecGetArrayF90(grid%xmol, x_p, ierr)
   endif
  
  ! call VecGetArrayF90(grid%vf, vf_p, ierr)   
   vf =0.D0
   if (grid%rk > 0.d0) then
      call VecGetArrayF90(grid%phis, phis_p, ierr)
    endif
!---------------------------------------------------------

 ! Open output file , write out title       
  if(grid%myrank==0)then 
    if (kplt < 10) then
      write(fname,'(a7,i1,a4)') 'pflow00', kplt, '.dat'
    else if (kplt < 100) then
      write(fname,'(a6,i2,a4)') 'pflow0', kplt, '.dat'
    else
      write(fname,'(a5,i3,a4)') 'pflow', kplt, '.dat'
    endif
    write(*,*) '--> write output file: ',fname
    open(unit=IUNIT3,file=fname,action="write")

    
    write(IUNIT3,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
    if( grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
                                   .or. grid%use_vadose == PETSC_TRUE &
                                   .or. grid%use_flash == PETSC_TRUE &
                                   .or. grid%use_richards == PETSC_TRUE) then
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'phase',q,'T',q,'p',q,'s(l)',q,'s(g)',q,&
       'U(l)',q,'U(g)',q,'xl(1)',q,'xl(2)',q,'xg(1)',q,'xg(2)',q,&
       'vf','"'
    else
 
   write(IUNIT3,'(''VARIABLES="'',3(a6,a3),a6,100(a3,a6))')"x",q,"y",q,"z",q,"ip",q,&
        'pressure',q,'temp', q,'sat',q,'conc',q,'vf','"' 
   endif
    write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz

   if  (grid%write_init == 1) then
         fname = 'pflow_init0.dat'
         write(*,*) '--> write output file: ',fname
         open(unit=1100,file=fname,action="write")
         if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE &
             .or. grid%use_flash == PETSC_TRUE & 
             .or. grid%use_vadose == PETSC_TRUE &
             .or. grid%use_richards == PETSC_TRUE) then
             write(1100,'(":kplt        steps     t       dt")')
             write(1100,'(2i10, 1p2e16.8)') kplt, grid%flowsteps, grid%t, grid%dt
             write(1100,'(":na  x   y   z          iphase           xx")') 
           else       
              write(1100,'(": i1  i2  j1  j2  k1  k2", &
                & "       p      ","      T      ","     sl(g)      ", &
                & "      C      ")')
          endif
   endif
  
  endif 

  !----------------------------------------------------------    
  !close(IUNIT3)
  !open(unit=IUNIT3,file=fname,action="write")
   if (grid%myrank == 0) &
   print *, 'pflow_output_new: tecplot format'         
   ndex=1
    do na=0, grid%nmax-1
      ifound=0
      if(grid%myrank==0) then
        ifound=0; vf=0.D0
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            jn= 1 + (n-1)*grid%ndof 
            nv = 1 + (n-1) * size_var_node
            if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
               .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
                iipha= iphase_p(n)! ; print *, iipha
                !print *,jn, xx_p(jn:jn+grid%ndof-1)
                xxx(:) = xx_p(jn:jn+grid%ndof-1)!; print *,jn, xxx
                vvar= var_p(nv:nv + size_var_use-1)
                if (grid%rk > 0.d0) vf = phis_p(n)
                !print *, na,n,iipha, xxx
             else
               pres(1:grid%nphase) = p_p(1+(n-1)*grid%nphase: n*grid%nphase) 
               temp = t_p(n)
               conc = c_p(n)
               sat (1:grid%nphase) = s_p(1+(n-1)*grid%nphase: n*grid%nphase)
               
!              print *,'pflow_output_new: ',grid%myrank,n,sat,temp,conc
               
               if (grid%use_2ph == PETSC_TRUE) xmol(1:grid%nphase) = x_p(1+(n-1)*grid%nphase: n*grid%nphase)
               if (grid%rk > 0.d0) vf = phis_p(n)
            endif  
            ndex=n
            exit
          endif
        enddo
      if(ifound==0)then  
        if  (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
               .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
            call MPI_Recv(iipha, 1, MPI_INTEGER, MPI_ANY_SOURCE, na+6553,PETSC_COMM_WORLD, status,ierr)  
            call MPI_Recv(xxx,grid%ndof, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
            na, PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(vvar, size_var_use ,MPI_DOUBLE_PRECISION, &
            MPI_ANY_SOURCE, na+na,PETSC_COMM_WORLD,status ,ierr)
            if (grid%rk > 0.d0) call MPI_Recv(vf,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
            na, PETSC_COMM_WORLD, status,ierr)
         else
            call MPI_Recv(pres,grid%nphase, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                 na, PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(temp,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                 na, PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(conc,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                 na, PETSC_COMM_WORLD, status,ierr)
            call MPI_Recv(sat,grid%nphase, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                 na, PETSC_COMM_WORLD, status,ierr)
            if (grid%use_2ph == PETSC_TRUE) &
               call MPI_Recv(xmol,grid%nphase, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
                 na, PETSC_COMM_WORLD, status,ierr)
           if (grid%rk > 0.d0)call MPI_Recv(vf,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
             na, PETSC_COMM_WORLD, status,ierr)
             
         endif      
         
       endif

        else ! other processors
          do n=ndex,grid%nlmax
            if(grid%nL2A(n) == na) then
              jn= 1 + (n-1)*grid%ndof 
              nv = 1 + (n-1) * size_var_node
             
             if  (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
               .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
               iipha=int(iphase_p(n))
              call MPI_Send(iipha, 1, MPI_INTEGER, 0, na+6553, PETSC_COMM_WORLD, ierr)   
              xxx= xx_p(jn:jn+grid%ndof-1)
              call MPI_Send(xxx, grid%ndof, MPI_DOUBLE_PRECISION,0, &
              na, PETSC_COMM_WORLD, ierr)  
              call MPI_Send(var_p(nv:nv + size_var_use-1), size_var_use, &
              MPI_DOUBLE_PRECISION, 0, na+na, PETSC_COMM_WORLD, ierr)   
              if (grid%rk > 0.d0)then
                 vf = phis_p(n)
                 call MPI_Send(vf, 1, MPI_DOUBLE_PRECISION,0, &
                               na, PETSC_COMM_WORLD, ierr)  
              endif
             else          
               pres(1:grid%nphase) = p_p(1+(n-1)*grid%nphase: n*grid%nphase) 
               temp = t_p(n)
               conc = c_p(n)
               sat (1:grid%nphase) = s_p(1+(n-1)*grid%nphase: n*grid%nphase)
               if (grid%use_2ph == PETSC_TRUE)&
                xmol(1:grid%nphase) = x_p(1+(n-1)*grid%nphase: n*grid%nphase)
               if (grid%rk > 0.d0) vf = phis_p(n) 
                
!             print *,'pflow_output_new-2: ',grid%myrank,n,sat,temp,conc
              
              call MPI_Send(pres, grid%nphase, MPI_DOUBLE_PRECISION,0, &
                   na, PETSC_COMM_WORLD, ierr)  
              call MPI_Send(temp, 1, MPI_DOUBLE_PRECISION,0, &
                   na, PETSC_COMM_WORLD, ierr)  
              call MPI_Send(pres, 1, MPI_DOUBLE_PRECISION,0, &
                   na, PETSC_COMM_WORLD, ierr)  
              call MPI_Send(sat, grid%nphase, MPI_DOUBLE_PRECISION,0, &
                   na, PETSC_COMM_WORLD, ierr)  
              if (grid%use_2ph == PETSC_TRUE) &
               call MPI_Send(xmol, grid%nphase, MPI_DOUBLE_PRECISION,0, &
                   na, PETSC_COMM_WORLD, ierr)  
              if (grid%rk > 0.d0)call MPI_Send(vf, 1, MPI_DOUBLE_PRECISION,0, &
              na, PETSC_COMM_WORLD, ierr)  
              endif                    
              ndex=n
              exit
            endif
          enddo  
        endif            

        call MPI_Barrier(PETSC_COMM_WORLD, ierr)
! End gather information


  
  if(grid%myrank ==0 )then
     if  (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
         .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
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
          vvar(3+ 7*grid%nphase: 2+ 7 *grid%nphase + grid%nphase* grid%nspec),& !Mol fractio
          vf
      else
        write(IUNIT3,'(1p10e12.4)') grid%x(na+1), grid%y(na+1), grid%z(na+1), &
          pres, temp, sat, conc, vf      
     endif
   


   if (grid%write_init == 1) then
              k= int(na/grid%nxy) + 1
              j= int(mod(na,grid%nxy)/grid%nx) + 1
              i= mod(mod(na,grid%nxy),grid%nx) + 1
     
           if (grid%use_2ph == PETSC_TRUE .or. &
                grid%use_mph == PETSC_TRUE .or. &
                grid%use_flash == PETSC_TRUE .or. &
                grid%use_vadose == PETSC_TRUE .or.&
                grid%use_richards == PETSC_TRUE ) then
                 write(1100,'(i10,1p10e16.8)')na,grid%x(na+1), grid%y(na+1), grid%z(na+1),&
                            real(iipha), xxx   
                 
            else
                  write(1100,'(6i4,1pe14.6,1p10e12.4)') i,i,j,j,k,k, &
                    pres, temp, sat, conc
            endif
   endif
   endif  ! endif of myrank     
 enddo  ! end loop of grid%nmax 

 
  if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE & 
        .or. grid%use_flash == PETSC_TRUE  .or. grid%use_richards == PETSC_TRUE)then
      call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
      call VecRestoreArrayF90(grid%xx, xx_p, ierr)
      call VecRestoreArrayF90(grid%var, var_p, ierr)
    else
      call VecRestoreArrayF90(grid%pressure, p_p, ierr)
      call VecRestoreArrayF90(grid%temp, t_p, ierr)
      call VecRestoreArrayF90(grid%conc, c_p, ierr)
      call VecRestoreArrayF90(grid%sat, s_p, ierr)
      if (grid%use_2ph == PETSC_TRUE) call VecRestoreArrayF90(grid%xmol, x_p, ierr)
   endif

 if (grid%rk > 0.d0) then
      call VecRestoreArrayF90(grid%phis, phis_p, ierr)
    endif

  if(grid%myrank ==0) then    
    close (IUNIT3)
    if(grid%write_init == 1) close(1100)
  endif 
 ! call MPI_BCast(kplt, 
 if (grid%iprint >= 1 .and. kplt > 0) then ! print out velocity field
   
    if(grid%myrank==0) then
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
  
   
      if (grid%nx > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlx0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlx', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=1001,file=fname,action="write")
        write(1001,'(''TITLE= "'',1pg12.4,'' [''a1,'']"'')') tyr,grid%tunit
        write(1001,'(''VARIABLES="'',a6,100(a3,a6))') &
            'x',q,'y',q,'z',q,'vlx','"'
        write(1001,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
   &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx-1,grid%ny,grid%nz
      endif
   
      if (grid%ny > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vly0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vly', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=1002,file=fname,action="write")
        write(1002,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(1002,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vly','"'
        write(1002,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny-1,grid%nz
       endif 
  
       if (grid%nz > 1) then
        if (kplt < 10) then
          write(fname,'(a10,i1,a4)') 'pflow_vlz0', kplt, '.dat'
        else
          write(fname,'(a9,i2,a4)') 'pflow_vlz', kplt, '.dat'
        endif
        write(*,*) '--> write output file: ',fname
        open(unit=1003,file=fname,action="write")
        write(1003,'(''TITLE= "'',1pg12.4,'' ['',a1,'']"'')') tyr,grid%tunit
        write(1003,'(''VARIABLES="'',a6,100(a3,a6))') &
          'x',q,'y',q,'z',q,'vlz','"'
        write(1003,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &        '' , J='',i4,'' , K='',i4)') tyr,grid%nx,grid%ny,grid%nz-1
       endif
     
  endif
   
   
   
!  print *,'Write out velocity field'
   call DAGlobalToLocalBegin(grid%da_3np_dof, grid%vl, INSERT_VALUES, &
                            vl_plot_loc, ierr)
   call DAGlobalToLocalend(grid%da_3np_dof, grid%vl, INSERT_VALUES, &
                            vl_plot_loc, ierr)
   
   call VecGetArrayF90(vl_plot_loc, vl_plot_loc_p,ierr)
   call VecGetArrayF90(grid%vl, vl_p,ierr)
!  print *,'velocity field got pointer' 
       ndex=1
  
      do k = 1, grid%nz
        do j = 1, grid%ny
          do i = 1, grid%nx
            na = i+(j-1)*grid%nx+(k-1)*grid%nxy -1
            vavel=0.D0
            velo=0.D0 
       !print *, grid%myrank,na, grid%nL2A(1), grid%nL2A(grid%nlmax)
      ifound=0
      if(grid%myrank==0) then
        ifound=0 
        if(na >= grid%nL2A(1) .and. na <= grid%nL2A(grid%nlmax))then
        do n=ndex,grid%nlmax
          if(grid%nL2A(n) == na) then
            ifound=1
            ng=grid%nL2G(n) -1
             
            
             if(grid%nx>1)then ! X-direction
            ip = 0; ng_upstream = ng - 1
            velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  
            if (i==1) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (i==grid%nx) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase) =.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)) 

            endif
           endif 

              if(grid%ny>1)then! y-direction
            ip = 1; ng_upstream = ng - grid%ngx
             velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  
            if (j==1 ) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (j==grid%ny ) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase) =.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)) 

            endif
           endif
           
            if(grid%nz>1)then  ! z-direction
            ip = 2; ng_upstream = ng - grid%ngxy
            velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  
            if (k==1) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (k==grid%nz ) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)=.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase))
          !   if(i== grid%nx) then
          !    print *,na,n,ng,ng_upstream, vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase),&
          !           vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase), vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase)
          !  endif          
     
            endif
           endif 
            ndex=n
            exit
          endif
        enddo
       endif
      if(ifound==0)then  
        call MPI_Recv(vavel,3*grid%nphase, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        call MPI_Recv(velo,3*grid%nphase, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        
       endif

    else
      if(na >= grid%nL2A(1) .and. na <= grid%nL2A(grid%nlmax))then
       do n=ndex,grid%nlmax
        if(grid%nL2A(n) == na) then
                         ifound=1
            ng=grid%nL2G(n) -1
            
             if(grid%nx>1)then
             ! X-direction
            ip = 0; ng_upstream = ng - 1
             velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  
            if (i==1) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (i==grid%nx) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)=.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)) 

            endif
         endif

            if(grid%ny>1)then            
               ! y-direction
            ip = 1; ng_upstream = ng - grid%ngx
            velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  
            if (j==1) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (j==grid%ny) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)=.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)) 

            endif
           endif

            if(grid%nz>1)then
             ! z-direction
            ip = 2; ng_upstream = ng - grid%ngxy
            velo(ip *grid%nphase +1:(ip +1)* grid%nphase) = &
               vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase: (ip+1)* grid%nphase + (n-1)*3*grid%nphase)  

            if (k==1) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase)*grid%tconv
             else if (k==grid%nz) then
              vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)&
              = vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)*grid%tconv
            else
             vavel(ip *grid%nphase +1:(ip +1)* grid%nphase)=.5D0 *grid%tconv * &
                (vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase:&
                    grid%nphase+ip*grid%nphase+(ng)*3*grid%nphase) +&
                 vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase:&
                 grid%nphase+ip*grid%nphase+(ng_upstream)*3*grid%nphase)) 
           ! if(i== grid%nx) then
           !   print *,na,n,ng,ng_upstream, vl_plot_loc_p(1+ip*grid%nphase+(ng)*3*grid%nphase),&
           !          vl_plot_loc_p(1+ip*grid%nphase+(ng_upstream)*3*grid%nphase),vl_p(1+ ip * grid%nphase + (n-1)*3*grid%nphase)
           ! endif          
            endif
           endif       
             call MPI_Send(vavel, 3*grid%nphase, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
             ndex=n
             call MPI_Send(velo, 3*grid%nphase, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)  
             ndex=n

          exit
        endif
      enddo
     endif 
    endif
                
    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
    
    if(grid%myrank ==0 )then
     write(IUNIT3,'(1p30e12.4)') grid%x(na+1),grid%y(na+1),grid%z(na+1), &
              (vavel(ip),ip=1,3*grid%nphase)

      if(grid%nx>1) then
         write(1001,'(1p10e12.4)') grid%x(na+1),grid%y(na+1),grid%z(na+1), &
              velo(1:grid%nphase) * grid%tconv
       endif 


      if(grid%ny>1) then
         write(1002,'(1p10e12.4)') grid%x(na+1),grid%y(na+1),grid%z(na+1), &
              velo(grid%nphase + 1: 2*grid%nphase)  *grid%tconv
       endif 
       
      if(grid%nz>1) then
         write(1003,'(1p10e12.4)') grid%x(na+1),grid%y(na+1),grid%z(na+1), &
              velo(2*grid%nphase+1: 3*grid%nphase) * grid%tconv
       endif 

    endif
  enddo
 enddo
enddo   
   call VecRestoreArrayF90(vl_plot_loc, vl_plot_loc_p,ierr)
       if(grid%myrank ==0 )then
          close(IUNIT3)
          if(grid%nx>1) close(1001)
          if(grid%ny>1) close(1002)
          if(grid%nz>1) close(1003)
       endif   
  endif   
 
  
   
     
 if (grid%iprint >= 2) call porperm_out (grid%t, grid%dt, grid%tconv, &
          kplt, grid%nx, grid%ny, grid%nz, grid%nmax, &
          grid%x, grid%y, grid%z, grid%flowsteps, grid%porosity,grid%perm_xx,&
          grid%myrank,grid%nlmax,grid%nL2A)
                        

 else    ! Switch to binary output using method in pflow_checkpoint
  !call DAGetCoordinate()     
  
  write(fname, '(a10,i6.6)') 'pflow.bin.', kplt
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, fname, FILE_MODE_WRITE, &
                             viewer, ierr)

  !--------------------------------------------------------------------
  ! Dump some important information such as simulation time, 
  ! time step size, etc.
  !--------------------------------------------------------------------
#if (PETSC_VERSION_RELEASE == 0 || PETSC_VERSION_SUBMINOR == 3)
  ! We manually specify the number of bytes required for the 
  ! checkpoint header, since sizeof() is not supported by some Fortran 
  ! compilers.  To be on the safe side, we assume an integer is 8 bytes.
  call PetscBagCreate(PETSC_COMM_WORLD, 88, bag, ierr)
  call PetscBagGetData(bag, header, ierr); CHKERRQ(ierr)

  ! Register variables that are passed into pflowGrid_step().
  call PetscBagRegisterInt(bag, header%ntstep, ntstep, "ntstep", &
                           "ntstep", ierr)
  call PetscBagRegisterInt(bag, header%kplt, kplt, "kplt", &
                           "kplt", ierr)
  call PetscBagRegisterInt(bag, header%iplot, iplot, "iplot", &
                           "iplot", ierr)
  call PetscBagRegisterInt(bag, header%iflgcut, iflgcut, "iflgcut", &
                           "iflgcut", ierr)
  call PetscBagRegisterInt(bag, header%ihalcnt, ihalcnt, "ihalcnt", &
                           "ihalcnt", ierr)
  call PetscBagRegisterInt(bag, header%its, its, "its", &
                           "its", ierr)
  
  ! Register relevant components of the pflowGrid.
  call PetscBagRegisterReal(bag, header%t, grid%t, "t", &
                            "Simulation time (years)", ierr)
  call PetscBagRegisterReal(bag, header%dt, grid%dt, "dt", &
                            "Current size of timestep (years)", ierr)
  call PetscBagRegisterInt(bag, header%flowsteps, grid%flowsteps, "flowsteps", &
                            "Total number of flow steps taken", ierr)
  call PetscBagRegisterInt(bag, header%kplot, grid%kplot, "kplot", &
                            "Printout steps", ierr)
  call PetscBagRegisterInt(bag, header%newtcum, grid%newtcum, "newtcum", &
                            "Total number of Newton steps taken", ierr)

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr)
  call PetscBagDestroy(bag, ierr)
#endif
  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  ! grid%xx is the vector into which all of the primary variables are 
  ! packed for the SNESSolve().
   
   
  ! here should dump out corrdinate, do not know how?
  

  ! If we are running with multiple phases, we need to dump the vector 
  ! that indicates what phases are present, as well as the 'var' vector 
  ! that holds variables derived from the primary ones via the translator.
  if(grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE .or. &
     grid%use_flash == PETSC_TRUE .or. grid%use_2ph == PETSC_TRUE .or. &
     grid%use_richards == PETSC_TRUE ) then
     call VecView(grid%xx, viewer, ierr)
     call VecView(grid%iphas, viewer, ierr)
    
     
    call VecView(var_plot, viewer, ierr)
     
    else
    call VecView(grid%pressure,viewer, ierr)
    call VecView(grid%temp,viewer, ierr)
    call VecView(grid%sat,viewer, ierr)
    call VecView(grid%conc,viewer, ierr)
    if(grid%use_2ph == PETSC_TRUE) call VecView(grid%xmol,viewer, ierr)
    call VecView(grid%hh, viewer, ierr)
    call VecView(grid%ddensity, viewer, ierr)
  endif  

  ! solid volume fraction
  if (grid%rk > 0.d0) then
    call VecView(grid%phis, viewer, ierr)
  endif

  ! Porosity and permeability.
  ! (We only write diagonal terms of the permeability tensor for now, 
  ! since we have yet to add the full-tensor formulation.)
  
if (grid%iprint >= 2) then  
  call VecView(grid%porosity, viewer, ierr)
  call VecView(grid%perm_xx, viewer, ierr)
  call VecView(grid%perm_yy, viewer, ierr)
  call VecView(grid%perm_zz, viewer, ierr)
endif

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr)

  if(grid%myrank == 0) write(*, '(a23, a16)') "Dumped binary file ", fname
    
 endif  ! endif of tecplot
   
        
 
    
  
  
  kplt = kplt + 1
  iplot = 0
 
 
 end subroutine pflow_output_new
 
   

  subroutine porperm_out(t, dt, tconv, kplt, nx, ny, nz, nmax, &
                         x, y, z, flowsteps, porosity, perm, myrank, nlmax,nla)
  
  use PetscRelWrappers  ! For petsc-release compatibility.

  implicit none


#include "definitions.h"
  
  integer, intent(in) :: kplt, nx, ny, nz, nmax, flowsteps, myrank
  real*8, intent(in) :: t, dt, tconv, x(:), y(:), z(:)
  Vec, intent(inout) :: porosity, perm
  
  integer :: ierr, nlmax,n, nla(:), status(MPI_STATUS_SIZE)
  
  
  real*8 :: fac, tyr, permm, poro
  
  real*8, pointer :: por_p(:), perm_p(:)
  
  character(len=20) :: fname
  character*3 :: q
  character*1 :: tab
  integer ifound, ndex, na


#include "definitions.h"

  tab = char(9)
  q = '","'
  tyr =t/tconv

! write out porosity, permeability

! porosity field
 
  if(myrank == 0) then
    
    if (kplt < 10) then
      write(fname,'(a10,i1,a4)') 'pflow_por0', kplt, '.dat'
    else
      write(fname,'(a9,i2,a4)') 'pflow_por', kplt, '.dat'
    endif
    
    write(*,*) '--> write output file: ',fname
    
    open(unit=IUNIT3,file=fname,action="write")
  endif
 
           
    call VecGetArrayF90(porosity, por_p, ierr)
    call VecGetArrayF90(perm, perm_p, ierr)
      
    fac = 1.d12
      

     if(myrank == 0)then
      write(IUNIT3,'(''TITLE= "'',1pg12.4,'' years"'')') tyr
      write(IUNIT3,'(''VARIABLES="'',a6,100(a3,a6))') &
      'x',q,'y',q,'z',q,'por',q,'perm','"' !,q,'permx',q,'permy',q,'permz','"'
      write(IUNIT3,'(''ZONE T= "'',1pg12.4,''",'','' I='',i4, &
 &    '' , J='',i4,'' , K='',i4)') tyr,nx,ny,nz
     endif
     
           
     ndex =1                  
    do na = 0, nmax-1
        
            ifound=0
      if(myrank==0) then
        ifound=0
        do n=ndex,nlmax
          if(nLA(n) == na) then
            ifound=1
            permm =perm_p(n)
            poro = por_p(n)
            ndex=n
            exit
          endif
        enddo
      if(ifound==0)then  
        call MPI_Recv(permm,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        call MPI_Recv(poro,1, MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, &
        na, PETSC_COMM_WORLD, status,ierr)
        endif

    else
      do n=ndex,nlmax
        if( nLA(n) == na) then
            permm =perm_p(n)
            poro = por_p(n)
          call MPI_Send(permm, 1, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)
          call MPI_Send(poro, 1, MPI_DOUBLE_PRECISION,0,na, PETSC_COMM_WORLD, ierr)    
          ndex=n
          exit
        endif
      enddo
    endif
                
    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
    
    if(myrank ==0 )then
        
        write(IUNIT3,'(1p10e12.4)') x(na+1), y(na+1), z(na+1), &
        poro,permm !, (log10(perm_p(i+3*(n-1))*fac),i=1,3)
     endif
    enddo
   
    call VecRestoreArrayF90(porosity, por_p, ierr)
    call VecRestoreArrayF90(perm, perm_p, ierr)
     if(myrank ==0 ) close (IUNIT3)

  end subroutine porperm_out

end module  pflow_output_module_new

    
       