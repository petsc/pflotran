 module translator_owg_module
 
  
 private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
  
   
! **3 phase condition*************************************************
! phase              Primary Variables      index
!  e                p, T, X(e,a), X(e,c)      1
!  g                p, T, X(g,a), X(g,c)      2 
!  l                p, T, X(l,a), X(l,c)      4
!  eg               p, T,  S(g),  X(g,c)      3 
!  el               p, T,  S(l),  X(l,c)      5
!  lg               p, T,  S(g),  X(g,c)      6
!  egl              p, T,  S(g),  S(l)        7
!**********************************************************************


! phase index 1.e; 2. l; 3. g
! within each phase component index : 1. H2O; 2. CO2; 3. Air

    
 public  pri_var_trans_owg_ninc,pri_var_trans_owg_winc , &
         translator_owg_check_phase_cond,translator_owg_step_maxchange,&
     translator_owg_massbal,translator_owg_switching
  
 real*8, private, parameter:: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                              fmwco2 = 44.0098D0, fmwoil= 142.D0
 real*8, private, parameter:: eps=5D-7, formeps = 5D-5
 real*8, private, parameter::yh2o_in_co2=1D-2   
 real*8, private :: tref
 real*8, private :: hkcoef(6)

 contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules


 subroutine translator_owg_massbal(grid)
 use pflow_gridtype_module
  implicit none
  type(pflowGrid) :: grid 
  
 
  integer :: ierr,icall
  integer :: n,n0,nc,np
  integer :: index, size_var_node
     
  PetscScalar, pointer ::  var_p(:),&
                           porosity_p(:), volume_p(:)
                           
!  integer, pointer ::iphase_p(:)
  
  real*8 ::  pvol,sum
  real*8, pointer ::  den(:),sat(:),xmol(:)
 
   real*8 :: tot(0:grid%nspec,0:grid%nphase), tot0(0:grid%nspec,0:grid%nphase)  
  data icall/0/

  call VecGetArrayF90(grid%var,var_p,ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
 
  size_var_node=(grid%ndof+1)*(2+7*grid%nphase +2*grid%nphase*grid%nspec)
  tot=0.D0

  do n = 1,grid%nlmax
    n0=(n-1)* grid%ndof
    index=(n-1)*size_var_node
    den=>var_p(index+3+grid%nphase: index+2+2*grid%nphase)
    sat=>var_p(index+2+1:index+2+grid%nphase)
    xmol=>var_p(index+2+7*grid%nphase+1:index+2+7*grid%nphase +&
         grid%nphase*grid%nspec)    
        
   
    pvol=volume_p(n)*porosity_p(n)
    
    do nc =1,grid%nspec
      do np=1,grid%nphase
        sum= sat(np)* xmol((np-1)*grid%nspec +nc)*den(np)
        tot(nc,np)= pvol*sum + tot(nc,np)
      !tot(0,np)=tot(0,np)+tot(nc,np)
      !tot(nc,0)=tot(nc,0)+tot(nc,np)
    enddo
    enddo
 nullify(sat, den,xmol) 
  enddo
 !  call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(grid%var,var_p,ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)

 
 
  if(grid%commsize >1)then
  do nc =0,grid%nspec
       do np=0,grid%nphase
    call MPI_REDUCE(tot(nc,np), tot0(nc,np),1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
!    call MPI_BCAST(tot0,(grid%nphase+1)*(grid%nspec+1),&
!        MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    enddo
  enddo
  if(grid%myrank==0) tot = tot0
   endif 
  
  
  if(grid%myrank==0)then
   do  np=1, grid%nphase
    tot(:,0)= tot(:,np) +  tot(:,0)
   enddo
 
    write(*,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
  & " gas:",1pe13.6, " tot:", 1p2e13.6," [kmol]")')&
  & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,3),tot(2,1)+tot(2,2)+tot(2,3)
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,*) 'time   dt   tot   totl   totg    toto'
      icall = 1
    endif
!   write(13,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
! &  " gas:",1pe13.6," tot:",1p2e13.6," [kmol]")')&
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2)
    write(13,'(1p100e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,&
 &  tot(2,:)
  endif   
  
  
  
  
 end subroutine translator_owg_massbal



 integer function translator_owg_check_phase_cond(iphase, var_node,num_phase,num_spec)
   implicit none
   integer iphase, num_phase, num_spec
   real*8, target:: var_node(:)
    
   integer ibase,succ,np,nc
   real*8, pointer :: t,p,satu(:),den(:), avgmw(:),h(:),u(:),pc(:),&
                      kvr(:),xmol(:),diff(:)
      
  real*8 sum     
    
  ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;         p=>var_node(ibase)
  ibase=ibase+1;         satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; avgmw=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; h=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; u=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; xmol=>var_node(ibase:ibase+num_phase*num_spec-1)
  ibase=ibase+num_phase*num_spec; diff=>var_node(ibase:ibase+num_phase*num_spec-1)
   
   succ=1
  if(iphase ==3)then
    do np =1, num_phase
    if(satu(np)>1D0 .or.satu(np)<0D0)then
     succ=-1  
       print *, 'phase=',iphase,satu(1:2)
     endif
    enddo
  endif
  
  if(iphase ==1 .or. iphase==2)then
    do np =1, num_phase
    sum=0D0
    do nc = 1, num_spec
     sum =sum + xmol((np-1)*num_spec+nc)
    enddo 
    if(sum > 1.D0+eps)then
     succ=-1 
     print *,'phase=',iphase,sum,xmol(1:4)
     endif 
    enddo
  endif
  nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)     
 translator_owg_check_phase_cond = succ
 end function translator_owg_check_phase_cond



 subroutine translator_owg_step_maxchange(grid)
   use pflow_gridtype_module
   type(pflowGrid), intent(inout) :: grid
  

  PetscScalar, pointer :: xx_p(:), yy_p(:), iphase_p(:),var_p(:),iphase_old_p(:)
  real*8 :: dsm,dcm, comp1,comp, cmp  
  real*8 :: dsm0,dcm0  
  integer n, j, iipha

   call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
   ! call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)
 
    grid%dtmpmax =0.D0
  
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)
  call VecGetArrayF90(grid%iphas_old, iphase_old_p,ierr)
  call VecGetArrayF90(grid%var, var_p, ierr)
  
  comp=0.D0;comp1=0.D0
  do  n=1, grid%nlmax
    n0=(n-1)*grid%ndof 
    iipha=int(iphase_p(n))
  if(abs(iipha-int(iphase_old_p(n)))<0.25D0 )then
    
    do j= 2, grid%ndof
    cmp=dabs(xx_p(n0+j)-yy_p(n0+j))
       
      select case(iipha)
      case(1,2,4)
          if(comp<cmp) comp=cmp
      case(3,5,6)
          if(j==2 .and. comp1<cmp)then
           comp1=cmp
        ! print *, '**max change s', comp1,cmp
        endif      
            if(j==3 .and. comp<cmp) comp=cmp
       case(7)    
        if(comp1<cmp) comp1=cmp
    case default
         print *, 'translator_owg_step_maxchange ::Error in phase assignment', iipha          
           stop
    end select    
     
       enddo
  
      else
    print *,'phase changed', n, iphase_p(n), iphase_old_p(n)
   endif
  enddo
  !call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  call VecRestoreArrayF90(grid%iphas_old, iphase_old_p,ierr)
  call VecRestoreArrayF90(grid%var, var_p, ierr)
 
  
 if(grid%commsize >1)then
    call MPI_ALLREDUCE(comp1, dsm0,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dsm0,1, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(comp, dcm0,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dcm0,1, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    comp1 = dsm0
    comp = dcm0
  endif 
  
  
  
        
   grid%dsmax=comp1
   grid%dcmax=comp
!   print *, 'max change',grid%dpmax,grid%dtmpmax,grid%dsmax,grid%dcmax
 end  subroutine translator_owg_step_maxchange




  subroutine Translator_OWG_Switching(xx,t,grid,icri,ichange,ierr)
  use pflow_gridtype_module
    use water_eos_module
    use gas_eos_module  
    use co2eos_module
    use span_wagner_module

  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  Vec, intent(in) :: xx
  integer icri,ichange, itable, ierr 

  PetscScalar, pointer :: xx_p(:), yy_p(:),iphase_p(:)
  integer*4 :: n,n0,index,ipr
  integer :: iipha,i 
  real*8 :: p2,p,tmp,t, xla
  real*8 :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  real*8 :: ug,  xphi,henry,co2_poyn
  real*8 :: xmol(grid%nphase*grid%nspec ),satu(grid%nphase) 
  real*8 :: x(1:grid%ndof)
  real*8 :: m11,m12, m21, m22, mb1, mb2, mm
  real*8 :: Henry_co2_oil, Henry_co2_water
  
 !print *, ' Translator_OWG_Switching begin'
! mphase code need assemble 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)
 
   
  ichange = 0   
  do n = 1,grid%nlmax
    ipr=0
    n0=(n-1)* grid%ndof
    iipha=iphase_p(n)
  x = xx_p(n0+1:n0+grid%ndof)
    p = x(1)  !  ; t= xx_p(n0+2)
 
 select case(iipha)
  
  case(1) ! only water phase
   p = x(1)
!   t = x(2)
   xmol(5)= x(2)
   xmol(6)= x(3)
   satu(1)=1.D0
   satu(2)=0.D0
   satu(3)=0.D0
   
   
  case(2) ! only Supercritical CO2
     p=x(1)
   xmol(5)= x(2)
   xmol(6)= x(3)   
    satu(1)=0.D0
   satu(2)=1.D0
   satu(3)=0.D0

     
  case(4) ! only  oil phase
    p=x(1)
   xmol(5)= x(2)
   xmol(6)= x(3)   
   satu(1)=0.D0
   satu(2)=0.D0
   satu(3)=1.D0
  
        
  
  case(3) ! water + SC 
  p=x(1)
  satu(1)=x(2)
  xmol(5)=x(3)
    satu(2)=1.D0-satu(1)
  satu(3)=0.D0
  
  
  case(5) ! water + oil
  p=x(1)
  satu(1)=x(2)
  xmol(5)=x(3)
  satu(2)=0.D0
  satu(3)=1.D0-satu(1)
  
  case(6) ! oil + SC  
  p=x(1)
  satu(2)=x(2)
  xmol(5)=x(3)
  satu(1)=0.D0
  satu(3)=1.D0-satu(2)
  
  case(7) ! 3 phase
  p=x(1)
  satu(1)=x(2)
  satu(2)=x(3)
  satu(3)=1.D0- satu(1)- satu(2)
   end select

 
!p2=p !*xmol(5)
 p2=p*xmol(5)
 if(p2>=5d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
     
      call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,grid%itable)
      dg= dg / fmwco2
      fg= fg * 1.D6 
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(p2,t,grid%scale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
      fg=p2
    endif
   xphi=fg/p2            

     call Henry_CO2_noderiv(xla,tmp,t,p*xmol(5),xphi,henry_co2_water,co2_poyn)

  
  hkcoef(1)= 100D0; hkcoef(4)=0.1D0; hkcoef(3)=0.01; hkcoef(6)=500D0
  hkcoef(2)= p/henry_co2_water
  hkcoef(5)= 1D1 * hkcoef(2)
          
     !xmol(2)= p*xmol(5)/ henry_co2_water 
     !xmol(8)= p*xmol(5)/ henry_co2_oil
   xmol(2) = xmol(5) * hkcoef(2)
   xmol(8) = xmol(5) * hkcoef(5)
   
!***********************************************************************
! Assumed simple solvablility
!  henry_oil_water = 0.0
    
! Or should call subroutines to determine oil solvability in water and SC phases here
!**************************************************************************

 !print *, 'finished xmol known'

   select case(iipha)
    case(1)    
   xmol(3) = xmol(6) *  hkcoef(3)
     xmol(1)=1.D0-xmol(2)-xmol(3)
   xmol(4) = xmol(1) /  hkcoef(1)
   xmol(9) = xmol(6)*  hkcoef(6)
   xmol(7) = xmol(4) * hkcoef(4)

   case(2) 
    xmol(4)= 1.D0- xmol(5)-xmol(6)  
    xmol(1)= xmol(4)*hkcoef(1)
    xmol(3)= xmol(6)*hkcoef(3)
    xmol(7)=xmol(4)*hkcoef(4)
    xmol(9)=xmol(6)*hkcoef(6)

   case(4) 
     xmol(9) = xmol(6) *hkcoef(6)
     xmol(7) = 1.D0- xmol(8) -xmol(9)
     xmol(4)= xmol(7)/ hkcoef(4) 
     xmol(6)= xmol(9)/hkcoef(6)
   xmol(1)= xmol(4) * hkcoef(1)
     xmol(3)= xmol(6) * hkcoef(3)
   
  case(3)  
    m11=hkcoef(1); m12= hkcoef(3); mb1=1.D0-xmol(2)
    m21=1D0; m22=1D0; mb2=1.D0-xmol(5)
    mm = m11*m22 -m12*m21  
      
    xmol(4)= (mb1*m22 -mb2*m12)/mm          
      xmol(6)= 1.D0- xmol(4)- xmol(5) 
    xmol(1) = xmol(4) *   hkcoef(1)
    xmol(3) = 1.D0 - xmol(1) -xmol(2)
    xmol(7) = xmol(4) * hkcoef(4)
    xmol(9) = xmol(6) * hkcoef(6)
    
     case(5)
      m11=hkcoef(1); m12= hkcoef(3); mb1=1.D0-xmol(2)
    m21=hkcoef(4); m22= hkcoef(6); mb2=1.D0-xmol(8)
    mm = m11*m22 -m12*m21  
      
    xmol(4)= (mb1*m22 -mb2*m12)/mm          
      xmol(6)=  (m11*mb2 - m21* mb1)/mm   
    xmol(1) = xmol(4) *   hkcoef(1)
    xmol(3) = 1.D0 - xmol(1) -xmol(2)
    xmol(7) = xmol(4) *   hkcoef(4)
    xmol(9) = 1.D0 -  xmol(7) -xmol(8)
    
   case(6)
      m11=1.D0; m12 =1.D0; mb1= 1.D0- xmol(5)
      m21=hkcoef(4); m22= hkcoef(6); mb2=1.D0-xmol(8)
    mm = m11*m22 -m12*m21 
    
    xmol(4)= (mb1*m22 -mb2*m12)/mm  
    xmol(6) = 1.D0 - xmol(4)- xmol(5)
    xmol(7)= hkcoef(4) *xmol(4)
      xmol(9) = 1.D0 -xmol(7) - xmol(8)
    xmol(1) = xmol(4) * hkcoef(1)
      xmol(3) = xmol(6) * hkcoef(3)
  ! print *, "/Case 6 ", n, x, xmol
    case(7)           
     m11= hkcoef(1)-hkcoef(2)
     m12= hkcoef(3)-hkcoef(2)
     mb1= 1.D0 - hkcoef(2)
     m21= hkcoef(4)-hkcoef(5)
     m22=  hkcoef(6)-hkcoef(5)
     mb2= 1.D0-hkcoef(5) 
     mm=m11*m22-m21*m12 
     
      xmol(4) =  (m22*mb1 - m12*mb2)/ mm
    xmol(6)=  ( m11*mb2 - m21 * mb1) /mm
    xmol(5)=1.D0 -xmol(4)- xmol(6)
    xmol(3)= hkcoef(3) * xmol(6)
      xmol(1)= xmol(4) * hkcoef(1)
    xmol(2)= 1.D0 - xmol(1)- xmol(3)
    xmol(7)=  xmol(4)* hkcoef(4)             
      xmol(9)= xmol(6)* hkcoef(6) 
    xmol(8) = 1.D0 - xmol(7)- xmol(9)
    
   end select
   
  ! print *, 'Switching: ',n, iipha, xmol 
   
 select case(icri)
  case(0)
   select case(iipha)
     case(2) !G
    if(xmol(4)* hkcoef(1)>=1.D0 .and. xmol(6)* hkcoef(6) >1.D0) then 
      write(*,'('' Gas -> Gas + W + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=formeps
      xx_p(n0+3)=1.D0- 2D0*formeps
      ichange = 1 ;ipr=1
     elseif(xmol(4)*hkcoef(1)>=1.D0)then
        write(*,'('' Gas -> Gas + W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 3
      xx_p(n0+2)= formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.D0)then 
        write(*,'('' Gas -> Gas + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 6
      xx_p(n0+2)=1.D0 - formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   case(1) !W
    if(xmol(5) >=1.025D0 .and. xmol(6)* hkcoef(6) >1.025D0) then 
         write(*,'('' W -> Gas + W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=1.D0- 2D0* formeps
      xx_p(n0+3)=formeps
      ichange = 1 ;ipr=1
     elseif(xmol(5)>=1.025D0)then
        write(*,'('' W -> Gas + W '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3),xmol
      iphase_p(n) = 3
      xx_p(n0+2)=1.D0 - formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.001D0)then 
        write(*,'('' W -> W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3),xmol
      iphase_p(n) = 5
      xx_p(n0+2)=1.D0  - formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   case(4) !O
    if(xmol(5) >=1.025D0 .and. xmol(4)* hkcoef(1) >1.001D0) then 
         write(*,'('' O -> Gas + W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=formeps
      xx_p(n0+3)=formeps
      ichange = 1 ;ipr=1
     elseif((xmol(5)+xmol(6))>=1.0D0)then
        write(*,'('' O -> Gas + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      iphase_p(n) = 6
      xx_p(n0+2)=formeps
      xx_p(n0+3)=0.85!(1.D0 - xmol(9)/hkcoef(6))
      ichange = 1 ;ipr=1
       elseif(xmol(4)*hkcoef(1) >1.0D0)then 
        write(*,'('' O -> W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      iphase_p(n) = 5
      xx_p(n0+2)=formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   case(3) !W+G
    if(satu(1)<=0.D0) then 
         write(*,'('' W + G -> G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 2
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(1)>1.D0)then
        write(*,'('' W + G -> W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 1
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.001D0)then 
     ! elseif((xmol(6)*hkcoef(6) + xmol(6)*hkcoef(6)+xmol(6)*hkcoef(6) )>1.001D0)then
      write(*,'('' W + G -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)= satu(1)-formeps
      xx_p(n0+3)= satu(2)-formeps
      ichange = 1 ;ipr=1
       endif
   case(5) ! W+ O
    if(satu(1)<=0.D0) then 
         write(*,'('' W + O -> O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 4
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(1)>1.D0)then
        write(*,'('' W + O -> W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 1
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
      elseif(xmol(5) >=.999D0)then 
   !   elseif((xmol(4)+xmol(5)+xmol(6)) >=1.025D0)then 
        write(*,'('' W + O -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=satu(1)
      xx_p(n0+3)= formeps
      ichange = 1 ;ipr=1
       endif
   case(6) !G + O
    if(satu(2)<=0.D0) then 
         write(*,'('' G + O -> O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 4
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(2)>=1.D0)then
        write(*,'('' G + O -> G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 2
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
       elseif(xmol(4)*hkcoef(1) >=1.001D0)then 
        write(*,'('' G + O -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)= formeps
      xx_p(n0+3)=satu(2)- formeps
      ichange = 1 ;ipr=1
       endif
    case(7) ! w+O
    if(satu(1)<=0.D0) then 
         write(*,'('' W + O +G -> O + G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 6
      xx_p(n0+2)=satu(2)
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       elseif(satu(3) <0.D0)then 
        write(*,'('' W + O +G -> W + G  '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 3
      xx_p(n0+2)=satu(1) 
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
     elseif(satu(2)<0.D0)then
        write(*,'('' W + O +G -> W +O'',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 5
      xx_p(n0+2)= satu(1)+ formeps
      xx_p(n0+3)= xmol(5)*.995
      ichange = 1 ;ipr=1
       endif
   end select

  case(1)
     select case(iipha)
     case(2) !G
     if(xmol(4)* hkcoef(1)>=1.D0 .and. xmol(6)* hkcoef(6) >1.D0) then 
      write(*,'('' Gas -> Gas + W + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=formeps
      xx_p(n0+3)=1.D0- 2D0*formeps
      ichange = 1 ;ipr=1
     elseif(xmol(4)*hkcoef(1)>=1.D0)then
        write(*,'('' Gas -> Gas + W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 3
      xx_p(n0+2)=formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.D0)then 
        write(*,'('' Gas -> Gas + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 6
      xx_p(n0+2)=1.D0 - formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   case(1) !W
    !if(xmol(5) >=1.0D0 .and. xmol(6)* hkcoef(6) >1.0D0) then 
       !  write(*,'('' W -> Gas + W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
     ! iphase_p(n) = 7
     ! xx_p(n0+2)=1.D0- 2D0* formeps
     ! xx_p(n0+3)=formeps
     ! ichange = 1 ;ipr=1
     !else
     if(xmol(5)>=1.D0)then
        write(*,'('' W -> Gas + W '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3),xmol
      iphase_p(n) = 3
      xx_p(n0+2)=1.D0 - formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.025D0)then 
        write(*,'('' W -> W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3),xmol
      iphase_p(n) = 5
      xx_p(n0+2)=1.D0  - formeps
      xx_p(n0+3)=xmol(5)*0.95
      ichange = 1 ;ipr=1
       endif
   case(4) !O
    if(xmol(5) >=1.0D0 .and. xmol(4)* hkcoef(1) >1.0D0) then 
         write(*,'('' O -> Gas + W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=formeps
      xx_p(n0+3)=formeps
      ichange = 1 ;ipr=1
     elseif((xmol(5)+xmol(6))>=1.0D0)then
        write(*,'('' O -> Gas + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      iphase_p(n) = 6
      xx_p(n0+2)=formeps
      xx_p(n0+3)=0.85!(1.D0 - xmol(9)/hkcoef(6))
      ichange = 1 ;ipr=1
       elseif(xmol(4)*hkcoef(4) >1.0D0)then 
        write(*,'('' O -> W + O '',i8,1p20e12.4)') n,xx_p(n0+1:n0+3), xmol
      iphase_p(n) = 5
      xx_p(n0+2)=formeps
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   case(3) !W+G
    if(satu(1)<=0.D0) then 
         write(*,'('' W + G -> G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 2
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(1)>1.D0)then
        write(*,'('' W + G -> W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 1
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
       elseif(xmol(6)*hkcoef(6) >1.0D0)then 
        write(*,'('' W + G -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)= satu(1)-formeps
      xx_p(n0+3)= satu(2)
      ichange = 1 ;ipr=1
       endif
   case(5) ! W+ O
    if(satu(1)<=0.D0) then 
         write(*,'('' W + O -> O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 4
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(1)>1.D0)then
        write(*,'('' W + O -> W '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 1
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
       elseif((xmol(4)+xmol(5)+xmol(6)) >=1.D0)then 
        write(*,'('' W + O -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)=satu(1)
      xx_p(n0+3)= formeps
      ichange = 1 ;ipr=1
       endif
   case(6) !G + O
    if(satu(2)<=0.D0) then 
         write(*,'('' G + O -> O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 4
      xx_p(n0+2)=xmol(5)
      xx_p(n0+3)=xmol(6)
      ichange = 1 ;ipr=1
     elseif(satu(2)>=1.D0)then
        write(*,'('' G + O -> G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 2
      xx_p(n0+2)= xmol(5)
      xx_p(n0+3)= xmol(6)
      ichange = 1 ;ipr=1
       elseif(xmol(4)*hkcoef(1) >=1.D0)then 
        write(*,'('' G + O -> W + G + O '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 7
      xx_p(n0+2)= formeps
      xx_p(n0+3)=satu(2)- formeps
      ichange = 1 ;ipr=1
       endif
    case(7) ! w+O
    if(satu(1)<=0.D0) then 
         write(*,'('' W + O +G -> O + G '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 6
      xx_p(n0+2)=satu(2)
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
     elseif(satu(2)<0.D0)then
        write(*,'('' W + O +G -> W +O'',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 5
      xx_p(n0+2)= satu(1)
      xx_p(n0+3)= xmol(5)
      ichange = 1 ;ipr=1
       elseif(satu(3) <0.D0)then 
        write(*,'('' W + O +G -> W + G  '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      iphase_p(n) = 3
      xx_p(n0+2)=satu(1) 
      xx_p(n0+3)=xmol(5)
      ichange = 1 ;ipr=1
       endif
   end select
  
    
 end select
 

 if(ipr ==1)then
!        iicap=int(icap_p(n))
!    iipha = int(iphase_p(n))
!    dif(1)= grid%difaq
!        dif(2)= grid%cdiff(int(ithrm_p(n)))
!        i=ithrm_p(n) 
   
!  call pri_var_trans_owg_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),iipha,&
 !       grid%scale,grid%nphase,grid%nspec,&
 !       iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
 !       grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
 !       grid%pcbetac(iicap),grid%pwrprm(iicap),dif,&
!    var_p((n-1)*size_var_node+1:(n-1)*size_var_node+size_var_use),&
!    grid%itable,ierr)

   
 !    index_var_begin=(n-1)*size_var_node+1
  !   index_var_end = index_var_begin -1 + size_var_use
     
   !  p1 = 1 + (n-1)*grid%ndof    
   !  call MPHASERes_ARCont(n, var_p(index_var_begin: index_var_end),&
!    porosity_p(n),volume_p(n),grid%dencpr(i), grid, Res, 0,ierr)

 !  print *,res,accum_p(p1:p1-1+grid%ndof)
   !accum_p(p1:p1-1+grid%ndof) =res
   ipr=0
   endif
  
  end do

  !print *,iphase_p
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
 ! print *,'Translator_OWG_Switching end'


  end subroutine Translator_OWG_Switching
  





  subroutine pri_var_trans_owg_ninc_3_3(x,tref,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
          var_node,itable,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    use water_eos_module
    use gas_eos_module  
    use co2eos_module
    use span_wagner_module
    use oil_eos_module
    use oil_pckr_module
  
    implicit none
    integer :: num_phase,num_spec,num_pricomp
    integer :: size_var_use
    real*8 energyscale, tref
    real*8,target:: var_node(:)
    integer :: iphase,itable,ierr
    integer :: ipckrtype !, ithrmtype
    
    real*8  :: pckr_sir(1:num_phase),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
    real*8  :: dif(1:num_phase)
  real*8 :: x(1:num_spec)
  real*8 :: m11,m12, m21, m22, mb1, mb2,mm

  
     
  real*8, pointer :: p,t
  real*8, pointer:: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
    real*8, pointer :: diff(:),xmol(:),satu(:)
  integer ibase 
    real*8 err
  
  real*8 p1,p2,tmp
  real*8 pw,dw_kg,dw_mol,hw,sat_pressure,vis_w,xphi
  real*8 dg,dddt,dddp,fg, dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  real*8 ug
  real*8 co2_phi, henry,co2_poyn
    real*8 stea,dsteamol,dstea_p,dstea_t, hstea,hstea_p,hstea_t,dstea
  real*8 den_oil, h_oil, visc_oil
  real*8 kr(num_phase), pckr_swir
  
  real*8 xla,vphi
    real*8 :: Henry_co2_oil, Henry_co2_water, x1, x2

  
  
   
   size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec 
    !print *, 'pri_var_trans_owg_3-3 begin', num_phase, num_spec, size_var_use
   
    ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;           p=>var_node(ibase)
  ibase=ibase+1;           satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; avgmw=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; h=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; u=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; xmol=>var_node(ibase:ibase+num_phase*num_spec-1)
  ibase=ibase+num_phase*num_spec; diff=>var_node(ibase:ibase+num_phase*num_spec-1)

  t=tref
   ! print *, 'pri_var_trans_owg_3-3 got pointer', iphase,xmol,satu,t,x
  select case(iphase)
  
  case(1) ! only water phase
   !print *, x,t
   p = x(1)
!   t = x(2)
   xmol(5)= x(2)
   xmol(6)= x(3)
   satu(1)=1.D0
   satu(2)=0.D0
   satu(3)=0.D0
   !print *, p,t, xmol,satu
   
  case(2) ! only Supercritical CO2
     p=x(1)
   xmol(5)= x(2)
   xmol(6)= x(3)   
    satu(1)=0.D0
   satu(2)=1.D0
   satu(3)=0.D0

     
  case(4) ! only  oil phase
     p=x(1)
   xmol(5)= x(2)
   xmol(6)= x(3)   
   satu(1)=0.D0
   satu(2)=0.D0
   satu(3)=1.D0
  
        
  
  case(3) ! water + SC 
  p=x(1)
  satu(1)=x(2)
  xmol(5)=x(3)
    satu(2)=1.D0-satu(1)
  satu(3)=0.D0
  
  
  case(5) ! water + oil
  p=x(1)
  satu(1)=x(2)
  xmol(5)=x(3)
  satu(2)=0.D0
  satu(3)=1.D0-satu(1)
  
  case(6) ! oil + SC  
  p=x(1)
  satu(2)=x(2)
  xmol(5)=x(3)
  satu(1)=0.D0
  satu(3)=1.D0-satu(2)
  
  case(7) ! 3 phase
  p=x(1)
  satu(1)=x(2)
  satu(2)=x(3)
  satu(3)=1.D0- satu(1)- satu(2)
   end select

  !print *, 'finished var '
 call oil_eos(t, p, 0.D0, 1.D0, den_oil, h_oil, energyscale, ierr)
 
!p2=p !*xmol(5)
 p2=p*xmol(5)
 if(p2>=5d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
     
      call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,itable)
      dg= dg / fmwco2
      fg= fg * 1.D6 
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
      fg=p2
    endif
   xphi=fg/p2            
     call Henry_CO2_noderiv(xla,tmp,t,p2,xphi,henry_co2_water,co2_poyn)

    ! henry_co2_oil = henry_co2_water * 1.D-1  
  
  hkcoef(1)= 100D0; hkcoef(4)=0.1D0; hkcoef(3)=0.01; hkcoef(6)=500D0
  hkcoef(2)= p/henry_co2_water
  hkcoef(5)= 1D1 * hkcoef(2)
          
     !xmol(2)= p*xmol(5)/ henry_co2_water 
     !xmol(8)= p*xmol(5)/ henry_co2_oil
   xmol(2) =xmol(5) * hkcoef(2)
   xmol(8)= xmol(5) * hkcoef(5)
   
!***********************************************************************
! Assumed simple solvablility
!  henry_oil_water = 0.0
    
! Or should call subroutines to determine oil solvability in water and SC phases here
!**************************************************************************

 !print *, 'finished xmol known'

   select case(iphase)
    case(1)    
   xmol(3) = xmol(6) *  hkcoef(3)
     xmol(1)=1.D0-xmol(2)-xmol(3)
   xmol(4) = xmol(1) /  hkcoef(1)
   xmol(9) = xmol(6)*  hkcoef(6)
   xmol(7) = xmol(4) * hkcoef(4)

   case(2) 
    xmol(4)= 1.D0- xmol(5)-xmol(6)  
    xmol(1)= xmol(4)*hkcoef(1)
    xmol(3)= xmol(6)*hkcoef(3)
    xmol(7)=xmol(4)*hkcoef(4)
    xmol(9)=xmol(6)*hkcoef(6)

   case(4) 
     xmol(9) = xmol(6) *hkcoef(6)
     xmol(7) = 1.D0- xmol(8) -xmol(9)
     xmol(4)= xmol(7)/ hkcoef(4) 
     xmol(6)= xmol(9)/hkcoef(6)
   xmol(1)= xmol(4) * hkcoef(1)
     xmol(3)= xmol(6) * hkcoef(3)

  case(3)  
    m11=hkcoef(1); m12= hkcoef(3); mb1=1.D0-xmol(2)
    m21=1D0; m22=1D0; mb2=1.D0-xmol(5)
    mm = m11*m22 -m12*m21  
      
    xmol(4)= (mb1*m22 -mb2*m12)/mm          
      xmol(6)= 1.D0- xmol(4)- xmol(5) 
    xmol(1) = xmol(4) *   hkcoef(1)
    xmol(3) = 1.D0 - xmol(1) -xmol(2)
    xmol(7) = xmol(4) * hkcoef(4)
    xmol(9) = xmol(6) * hkcoef(6)
    
     case(5)
      m11=hkcoef(1); m12= hkcoef(3); mb1=1.D0-xmol(2)
    m21=hkcoef(4); m22= hkcoef(6); mb2=1.D0-xmol(8)
    mm = m11*m22 -m12*m21  
      
    xmol(4)= (mb1*m22 -mb2*m12)/mm          
      xmol(6)=  (m11*mb2 - m21* mb1)/mm   
    xmol(1) = xmol(4) *   hkcoef(1)
    xmol(3) = 1.D0 - xmol(1) -xmol(2)
    xmol(7) = xmol(4) *   hkcoef(4)
    xmol(9) = 1.D0 -  xmol(7) -xmol(8)
    
   case(6)
      m11=1.D0; m12 =1.D0; mb1= 1.D0- xmol(5)
      m21=hkcoef(4); m22= hkcoef(6); mb2=1.D0-xmol(8)
    mm = m11*m22 -m12*m21 
    
    xmol(4)= (mb1*m22 -mb2*m12)/mm  
    xmol(6) = 1.D0 - xmol(4)- xmol(5)
    xmol(7)= hkcoef(4) *xmol(4)
      xmol(9) = 1.D0 -xmol(7) - xmol(8)
    xmol(1) = xmol(4) * hkcoef(1)
      xmol(3) = xmol(6) * hkcoef(3)
  ! print *, 'case 6: ', x, xmol
    case(7)           
     m11= hkcoef(1)-hkcoef(2)
     m12= hkcoef(3)-hkcoef(2)
     mb1= 1.D0 - hkcoef(2)
     m21= hkcoef(4)-hkcoef(5)
     m22=  hkcoef(6)-hkcoef(5)
     mb2= 1.D0-hkcoef(5) 
     mm=m11*m22-m21*m12 
     
      xmol(4) =  (m22*mb1 - m12*mb2)/ mm
    xmol(6)=  ( m11*mb2 - m21 * mb1) /mm
    xmol(5)=1.D0 -xmol(4)- xmol(6)
    xmol(3)= hkcoef(3) * xmol(6)
      xmol(1)= xmol(4) * hkcoef(1)
    xmol(2)= 1.D0 - xmol(1)- xmol(3)
    xmol(7)=  xmol(4)* hkcoef(4)             
      xmol(9)= xmol(6)* hkcoef(6) 
    xmol(8) = 1.D0 - xmol(7)- xmol(9)
      !  print *, 'w+G+O xmol:',xmol
   end select
! print *, 'finished xmol assign', xmol
 
! Then for the 
    avgmw(1)= xmol(1)* fmwh2o + xmol(2) * fmwco2 + xmol(3) * fmwoil
    avgmw(2)= xmol(4)* fmwh2o + xmol(5) * fmwco2 + xmol(6) * fmwoil         
  avgmw(3)= xmol(7)* fmwh2o + xmol(8) * fmwco2 + xmol(9) * fmwoil  
  
  diff(1:num_spec)=dif(1)
    diff(num_spec +1: 2*num_spec)=dif(2)  
    diff(2*num_spec +1: 3*num_spec)=dif(3)        
               
! Pc---kr coorelation                  

   call oil_pckr_noderiv(ipckrtype,pckr_sir(1),pckr_sir(3),pckr_lambda,pckr_alpha,&
        pckr_m,pckr_pcmax,satu(1),satu(3),pc,kr,pckr_betac,pckr_pwr)



! Water phase **********************************     
  pw=p  
    call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)                 
     call PSAT(t, sat_pressure, ierr)
  call VISW_noderiv(t,pw,sat_pressure,vis_W,ierr)
    call oil_eos(t,p, 0.D0, 1.D0, den_oil, h_oil, energyscale, ierr)
  call Vis_oil(p,t,visc_oil,ierr)
  den_oil=dw_kg/fmwoil


 ! Garcia mixing
  ! if((xmol(1)+ xmol(2))>1D-4)then
   if(iphase/=4 .and. iphase/=6)then
   x1=xmol(1); x2=xmol(2)
   tmp =  x1* fmwh2o + x2 * fmwco2
   vphi=1D-6*(37.51D0+t*(-9.585D-2+t*(8.74D-4-t*5.044D-7)))
   den(1)=dw_mol*fmwh2o/(1D0-(fmwco2*1D-3-dw_mol*fmwh2o*vphi)*x2/(tmp*1D-3))
   den(1)=den(1)/tmp
  ! den(1)=1.D0/((xmol(1)+xmol(2))/den(1)+xmol(3)/den_oil)
   else
   den(1)=dw_mol
   endif
   
!  den(1)=  1.D0/( xmol(1)/dw_mol+ xmol(2)/dg + xmol(3)/den_oil)
    h(1) = hw * xmol(1) + hg*xmol(2) + h_oil*xmol(3) 
    u(1) = h(1) - pw /den(1)* energyscale 
  kvr(1)=kr(1)/vis_w
 
 
 ! oil phase **************************************
   ! Ideal mixing for density
    tmp = vis_w * xmol(7) + visg*xmol(8) + visc_oil*xmol(9)
  den(3)=1.D0/(xmol(7)/dw_mol + xmol(8)/dg + xmol(9)/den_oil)
    h(3)= hw * xmol(7) + hg*xmol(8) + h_oil*xmol(9)
  u(3) =  h(3) - pw /den(3)* energyscale
  kvr(3)=kr(3)/tmp

! SC phase *******************************************
    tmp = vis_w * xmol(4) + visg*xmol(5) + visc_oil*xmol(6)
  den(2)= 1.D0/(xmol(4)/dw_mol + xmol(5)/dg + xmol(6)/den_oil)
    h(2)= hw * xmol(4) + hg*xmol(5) + h_oil*xmol(6)
    u(2) =  h(2) - pw /den(2)* energyscale
    kvr(2)=kr(2)/ tmp
  !print *, 'trans vis', kvr, xmol
  
  nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
 end subroutine pri_var_trans_owg_ninc_3_3
  



! **2 phase condition**************************************************
! phase             Primary Variables      index
!   e                p, T, X(e,c)                  1
!   g                p, T, X(g,a)                  2 
!   eg               p, T, S(g)                    3
!**********************************************************************



 
 subroutine pri_var_trans_owg_ninc(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
          var_node,itable,ierr,phi_co2, tref)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    implicit none
    integer :: num_phase,num_spec,num_pricomp
    integer :: size_var_use
  real*8 x(:),energyscale
    real*8 var_node(1:2 + 7*num_phase + 2* num_phase*num_spec)
  real*8 :: dif(:)
  integer ::iphase, itable,ierr
  integer :: ipckrtype !, ithrmtype
     
    
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr 
    real*8 :: phi_co2
  
  real*8 :: xphi_co2=1.D0
  real*8 :: tref
  

    size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
    if((num_phase == 3).and.( num_spec ==3)) then
   ! print *, 'pri_var_trans_owg_ninc begin 3-3',x, tref,iphase
     call pri_var_trans_owg_ninc_3_3( x,tref,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
          var_node,itable,ierr)
  !print *, 'pri_var_trans_owg_ninc end 3-3',  var_node      
    else 
   print *, 'Wrong phase-specise combination. Stop.'
   stop
   endif
  end subroutine pri_var_trans_owg_ninc   
  
  
 subroutine pri_var_trans_owg_winc(x,delx,iphase,energyscale,num_phase,num_spec,&
                    num_dof, ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
          var_node,itable,ierr, tref)
  integer :: num_phase,num_spec, num_dof
  integer :: size_var_use,size_var_node
    

    real*8 x(:),delx(:),energyscale
    real*8 var_node(:)
  real*8 :: dif(:)
  integer ::iphase,itable,ierr
  integer :: ipckrtype !, ithrmtype
    real*8  tref
   
    

    real*8 xx(1:num_spec), tmp
    real*8 pckr_sir(1: num_phase),pckr_lambda,pckr_alpha,pckr_m,&
           pckr_pcmax,pckr_betac,pckr_pwr

    size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
    size_var_node = (num_spec+1)*size_var_use
  !print *, 'pri_var_trans_owg_winc Begin' 
  if(num_dof /= num_spec) then
     print *, "Wrong combination:: STOP!!!"
     stop
  endif   
   do n=1, num_dof
         xx=x;  xx(n)=x(n)+ delx(n)
  ! note: var_node here starts from 1 to grid%ndof*size_var_use
     ! print *, 'pri_var_trans_owg_winc ', n, x,xx 
      call pri_var_trans_owg_ninc(xx,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
          var_node((n-1)*size_var_use+1:n*size_var_use),itable,ierr, tmp,tref)
  
      enddo                
  !print *, 'pri_var_trans_owg_winc End'
          
 end subroutine pri_var_trans_owg_winc
  
end module translator_owg_module
