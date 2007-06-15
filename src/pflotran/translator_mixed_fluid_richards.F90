module translator_Richards_module
 
  
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
!#ifdef USE_PETSC216
!#include "include/finclude/petscsles.h"
!#endif
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
  

! **3 phase condition*************************************************
! phase                             Primary Variables      index
!  e                p, T, X(e,a), X(e,c)          1
!  g                               p, T, X(g,a), X(g,c)          2 
!   l                p, T, X(l,a), X(l,c)      4  
!   eg                p, T,  S(g),  X(g,c)        3 
!   el                p, T,  S(l),  X(l,c)          5
!   lg                p, T,  S(g),  X(g,c)          6
!   egl                p, T,  S(g),  S(l)            7
!**********************************************************************


! phase index 1.e; 2. l; 3. g
! within each phase component index : 1. H2O; 2. CO2; 3. Air

    
  public  pri_var_trans_Richards_ninc,pri_var_trans_Richards_winc ,translator_Ric_step_maxchange, &
          translator_Richards_get_output,translator_check_cond_Richards,translator_Richards_massbal, &
          Translator_Richards_Switching
     
     
  real*8, private, parameter :: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                                fmwco2 = 44.0098D0
  real*8, private, parameter :: eps=5D-7 , formeps=5D-5
  real*8, private, parameter :: rgasj   = 8.3143    ![J/K/mol]

contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules


subroutine translator_Richards_massbal(grid)
 
  use pflow_gridtype_module
  
  implicit none

  type(pflowGrid) :: grid 
  
 
  integer :: ierr,icall
  integer :: n,n0,nc,np,n2p,n2p0
  real*8 x,y,z,nzm,nzm0, nxc,nxc0,c0, c00,nyc,nyc0,nzc,nzc0,nsm,nsm0,sm 
  integer :: index, size_var_node
     
  PetscScalar, pointer ::  var_p(:),&
                           porosity_p(:), volume_p(:)
                           
  PetscScalar, pointer ::iphase_p(:)
  
  real*8 ::  pvol,sum
  real*8, pointer ::  den(:),sat(:),xmol(:)
 
  real*8 :: tot(0:grid%nspec,0:grid%nphase), tot0(0:grid%nspec,0:grid%nphase)  
  data icall/0/

  call VecGetArrayF90(grid%var,var_p,ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
 
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
! print *,nzc,c0
  enddo
 !  call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(grid%var,var_p,ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
 
 
  if (grid%commsize >1) then
          
    do nc = 0,grid%nspec
      do np = 0,grid%nphase
        call MPI_REDUCE(tot(nc,np),tot0(nc,np),1,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,0,PETSC_COMM_WORLD,ierr)
   
!       call MPI_BCAST(tot0,(grid%nphase+1)*(grid%nspec+1),&
!            MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
      enddo
    enddo
  endif 
   if (grid%myrank==0) tot = tot0

  
  if (grid%myrank==0) then
  
  
    write(*,'(" Total H2O: liq:",1p, e13.6," gas:",1p, e13.6, " tot:", 1p, &
  &           2e13.6, " [kmol]",1p, 3e13.6)') &
          tot !,nzc,nzm,nsm
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,*) '# time   dt   totl   totg   tot n2p'
      icall = 1
    endif
!   write(13,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
! &  " gas:",1pe13.6," tot:",1p2e13.6," [kmol]")')&
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2)
    write(13,'(1p9e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,tot(:,1)
  endif    
  
end subroutine translator_Richards_massbal


integer function translator_check_cond_Richards(iphase, var_node,num_phase,num_spec)

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
  if (iphase == 3) then
    do np =1, num_phase
      if (satu(np)>1D0 .or.satu(np)<0D0) then
        succ = -1  
        print *, 'phase=',iphase,satu
      endif
    enddo
  endif
  
   nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
 translator_check_cond_Richards = succ

end function translator_check_cond_Richards


subroutine translator_Richards_get_output(grid)

  use pflow_gridtype_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  integer :: ierr
  
  PetscScalar, pointer :: t_p(:),p_p(:),c_p(:),s_p(:),cc_p(:),var_P(:)
  integer :: n, index_var_begin ,jn, size_var_node
  PetscScalar, pointer :: p,t,satu(:),xmol(:)
  
  call VecGetArrayF90(grid%var, var_p, ierr)
  call VecGetArrayF90(grid%pressure, p_p, ierr)
  call VecGetArrayF90(grid%temp, t_p, ierr)
  call VecGetArrayF90(grid%xmol, c_p, ierr)
  call VecGetArrayF90(grid%sat, s_p, ierr)
  call VecGetArrayF90(grid%conc, cc_p, ierr)
  !print *,' translator_mph_get_output gotten pointers'
  
  size_var_node=(grid%ndof+1)*(2+7*grid%nphase +2*grid%nphase*grid%nspec)
  
  do n = 1, grid%nlmax
    index_var_begin = (n-1) * size_var_node
    jn = 1 + (n-1)*grid%nphase 
    
    p_p(jn) = var_p(index_var_begin + 2)! - var_p(index_var_begin+5*grid%nphase+3)
   ! p_p(jn+1) = var_p(index_var_begin + 2) - var_p(index_var_begin+5*grid%nphase+4)
   
    t_p(n) = var_p(index_var_begin + 1)

    c_p(jn) = 0.D0!var_p(index_var_begin+7*grid%nphase+4)
    if(grid%nspec>1)  c_p(jn) = var_p(index_var_begin +2+ 7*grid%nphase +2)
   ! c_p(jn+1) = var_p(index_var_begin+7*grid%nphase+6)
 !   cc_p(n) = c_p(jn+1)
  
    s_p(jn) = var_p(index_var_begin + 3) 
 !   s_p(jn+1)=1.D0 -  s_p(jn)
  enddo
 
  call VecRestoreArrayF90(grid%var, var_p, ierr)
  call VecRestoreArrayF90(grid%pressure, p_p, ierr)
  call VecRestoreArrayF90(grid%temp, t_p, ierr)
  call VecRestoreArrayF90(grid%xmol, c_p, ierr)
  call VecRestoreArrayF90(grid%sat, s_p, ierr)
  call VecRestoreArrayF90(grid%conc, cc_p, ierr)
 
! work only for 2 phases
end subroutine translator_Richards_get_output


subroutine translator_Ric_step_maxchange(grid)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  

  PetscScalar, pointer :: xx_p(:), yy_p(:), iphase_p(:),var_p(:),iphase_old_p(:)
  real*8 :: dsm,dcm, comp1,comp, cmp  
  real*8 :: dsm0,dcm0  
  integer :: n, j, n0, ierr

  call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
  call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
  call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)

   grid%dsmax=0.D0
   grid%dcmax=0.D0
   if (grid%myrank == 0) &
     print *, 'ric max change',grid%dpmax,grid%dtmpmax,grid%dsmax,grid%dcmax
end subroutine translator_Ric_step_maxchange


subroutine Translator_Richards_Switching(xx,grid,icri,ichange)

  use pflow_gridtype_module
  use water_eos_module
  use gas_eos_module  
    
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  Vec, intent(in) :: xx
  integer :: icri,ichange 

  PetscScalar, pointer :: xx_p(:), yy_p(:),iphase_p(:)
  integer :: n,n0,index,ipr
  integer :: ierr,iipha,i 
  real*8 :: p2,p,tmp,t, xla, sat_pressure
  real*8 :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  real*8 :: ug,xphi,henry,co2_poyn
  real*8 :: xmol(grid%nphase*grid%nspec ),satu(grid%nphase) 

! mphase code need assemble 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)

  ichange = 0
  do n = 1,grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ipr=0
    n0=(n-1)* grid%ndof
    iipha=iphase_p(n)
    p = xx_p(n0+1); t= xx_p(n0+2)
    
     select case(iipha) 
      case(1) 
        xmol(2)= xx_p(n0+3)
        xmol(1)=1.D0 - xmol(2)
        satu(1)=1.D0; satu(2)=0.D0
      case(2)  
        xmol(4)= xx_p(n0+3)
        xmol(3)=1.D0 - xmol(4)
        satu(1)=0.D0; satu(2)=1.D0
      case(3) 
        satu(2)= xx_p(n0+3) 
        satu(1)= 1.D0- satu(2)
        xmol(3)= sat_pressure /p ; xmol(4)=1.D0-xmol(3)
    end select
  
    p2=p*xmol(4)
!    p2=p*xmol(4)
    
    call ideal_gaseos_noderiv(p ,t,grid%scale,dg,hg,ug)
    call visgas_noderiv(t,p2,p,dg,visg)
    fg = p2
   
    xphi = 1.D0
!   xphi = 1.d0 ! test-pcl

!   call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
!  call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
    call   Henry_air_noderiv(p,t,sat_pressure , Henry)   
   !  henry= henry/p 
    !note: henry = H/phi
    
!    print *,'translator_mixed: ',iipha,xla,tmp,t,p,xmol(4),p*xmol(4), &
!    xphi,henry,dg,fg,hg
     
  !   tmp = (p-sat_pressure)/(henry-sat_pressure)
 !    err= dabs(tmp-xmol(4))
!   tmp=xmol(4)
 ! enddo
  
      
 !print *, 'out 2 phase solver'
    select case(iipha)     
      case(1)
        xmol(4)=xmol(2)*henry/p   
    ! if(iphase /= 1)then
        xmol(3)= sat_pressure *xmol(1)/p
  !if(xmol(3)<0.D0) xmol(3)=0.D0
      case(2)   
        xmol(2)= p*xmol(4)/henry
    !xmol(2)= 1.D0/(1.D0 + 1.D0 /(henry * xmol(4) * p* 1D-5 *xphi * grid%fmwh2o ))
        xmol(1)=1.D0-xmol(2)
        if (xmol(1)<0.D0) xmol(1)=0.D0
      case(3)
        xmol(1)=(Henry - p)/(Henry -sat_pressure)
        xmol(2)= 1.D0- xmol(1)
        xmol(3)= xmol(1) * sat_pressure / p
        xmol(4)= 1.D0-xmol(3)
    end select

    select case(icri)
      case(0)
        if (xmol(3) > sat_pressure/p .and. iipha==2 ) then
          write(*,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 3
          xx_p(n0+3)=1D0-eps
          ichange = 1 ;ipr=1
!       xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
        endif
       ! gas ->  2ph 

    !if (xx_p(n0+3) > 1.025D0  .and. iipha==1) then
        if ((xmol(4)+ xmol(3)) > 1.025D0  .and. iipha==1) then
          write(*,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3), &
                                                     xmol(3)+xmol(4),xmol(3), &
                                                     xmol(4)
      !write(IUNIT2,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 3
     
     !tmp= (xmol(4)-1.D0)/(xmol(4)/xmol(2))*den(1)/den(2)
    !if(tmp>eps) tmp=eps
     ! xx_p(n0+3)=max(tmp,formeps)
          xx_p(n0+3)=formeps
!     xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
    !if(dabs(xmol(4)-1.D0) < eps) xx_p(n0+1)=xx_p(n0+1)* (1.D0-eps)
          ichange = 1;ipr=1
        endif
      
        if (satu(2)>1.0D0.and. iipha==3 ) then
  !if(xx_p(n0+3)> 1.D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
       ! write(IUNIT2,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 2
          xx_p(n0 + 3) = xmol(4)
!    xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)  
          ichange =1    ;ipr=1  
        endif
    

   ! if(sat(2)<= -formeps .and. iipha==3 )then
        if (satu(2)<= 0.0D0 .and. iipha==3 ) then
          write(*,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3), &
                                                     satu(1),satu(2)
      ! write(IUNIT2,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 1 ! 2ph -> Liq
          ichange = 1;ipr=1
          tmp = xmol(2) * 0.9995
          xx_p(n0 + 3)=tmp
!    xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
    !xx_p(n0+2) =  xx_p(n0+2)*(1.D0-eps)
        endif
     
    
!  if(ichange ==1) then
!    call 
      case(1)
  
        if ( xmol(3)> sat_pressure/p .and. iipha==2 ) then
          write(*,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 3
          xx_p(n0+3)=1D0-eps
          ichange = 1 ;ipr=1
        !  xx_p(n0+2)= yy_p(n0+2)
        endif
       ! gas ->  2ph 

        if ( (xmol(4)+xmol(3)) > 1.0D0  .and. iipha==1) then
          write(*,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3), &
                                                     xmol(3)+xmol(4),xmol(3), &
                                                     xmol(4)
      !write(IUNIT2,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 3
     
     !tmp= (xmol(4)-1.D0)/(xmol(4)/xmol(2))*den(1)/den(2)
    !if(tmp>eps) tmp=eps
     ! xx_p(n0+3)=max(tmp,formeps)
          xx_p(n0+3)=formeps
    !if(dabs(xmol(4)-1.D0) < eps) xx_p(n0+1)=xx_p(n0+1)* (1.D0-eps)
          ichange = 1;ipr=1
        endif
      
        if (satu(2)>1.0D0.and. iipha==3 ) then
  !if(xx_p(n0+3)> 1.D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
       ! write(IUNIT2,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 2
          xx_p(n0 + 3) = xmol(4)  
          ichange =1; ipr=1  
        endif
    

   ! if(sat(2)<= -formeps .and. iipha==3 )then
        if (satu(2)<= 0.0D0 .and. iipha==3 ) then
          write(*,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3), &
                                                     satu(1),satu(2)
      !write(IUNIT2,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_p(n) = 1 ! 2ph -> Liq
          ichange = 1;ipr=1
          tmp = xmol(4) * 0.9995
          xx_p(n0 + 3)=tmp
    !xx_p(n0+2) =  xx_p(n0+2)*(1.D0-eps)
        endif

    end select


    if (ipr ==1) then
!        iicap=int(icap_p(n))
!    iipha = int(iphase_p(n))
!    dif(1)= grid%difaq
!        dif(2)= grid%cdiff(int(ithrm_p(n)))
!        i=ithrm_p(n) 

!  call pri_var_trans_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),iipha,&
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

end subroutine Translator_Richards_Switching
  
! **richards phase condition**************************************************
! phase                             Primary Variables      index
!   e                p, T, X(e,c)                  1
!   g                p, T, X(g,a)                  2 
!   eg                              p, T, S(g)                    3
!**********************************************************************
subroutine pri_var_trans_Richards_ninc_2_2(x,iphase,energyscale,num_phase,num_spec,&
                                      ipckrtype,pckr_sir,pckr_lambda, &
                                      pckr_alpha,pckr_m,pckr_pcmax,pckr_betac, &
                                      pckr_pwr,dif,var_node,itable,ierr, pref)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  use water_eos_module
  use gas_eos_module  
  use pckr_module
  
  implicit none

  integer :: num_phase,num_spec, itable, ierr
  integer :: size_var_use 
  real*8 :: x(1:num_spec+1),energyscale
  real*8, target:: var_node(:)
  integer ::iphase
  integer :: ipckrtype !, ithrmtype
  real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
  real*8 :: dif(:)

   
 !   integer size_var_node = (grid%ndof+1)*size_var_use

  real*8, pointer :: t ,p
  real*8, pointer :: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
  real*8, pointer :: xmol(:),satu(:),diff(:)
  integer :: ibase 
  
  real*8 :: p1,p2,tmp, pref
  real*8 :: pw,dw_kg,dw_mol,hw,sat_pressure,visl,xphi, dco2
  real*8 :: dg,dddt,dddp,fg, dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  real*8 :: ug
  real*8 :: co2_phi, henry,co2_poyn
  real*8 :: stea,dsteamol,dstea_p,dstea_t, hstea,hstea_p,hstea_t,dstea
  real*8 :: kr(num_phase), pckr_swir
  real*8 :: err,xla,vphi

  
  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  !pckr_swir=pckr_sir(1)
  
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

!note: xmol means molar concentration, room for multiple tracers test. 

 satu =0.D0
 h=0.D0
 u=0.D0
 den=0.D0
 avgmw =0.D0
 xmol =0.D0
 kr=0.D0
 diff =0.D0
 
 p=x(1)  
 t=x(2)
 
 pc(1) = pref - p
 xmol(1)=1.D0
 if(num_spec>1) xmol(2:num_spec)=x(3: num_spec +1)   

!***************  Liquid phase properties **************************
  avgmw(1)=  fmwh2o 

 ! no more calculation on gas phase
!  avgmw(2)= xmol(3)* fmwh2o + xmol(4) * fmwa 
  pw=pref   
   if(pc(1)>0.D0)then
    iphase = 3
    call pflow_pckr_richards(ipckrtype,pckr_sir(1),pckr_lambda,pckr_alpha, &
                            pckr_m,pckr_pcmax,satu(1),pc,kr,pckr_betac, &
                            pckr_pwr)
  else
    iphase = 1
    pc(1)=0.D0
    satu(1)= 1.D0  
    kr(1)=1.D0    
    pw=p
  endif  

  call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
   !    call VISW(t,pw,sat_pressure,visl,tmp,tmp2,ierr)
   ! call VISW_FLO(t,dw_mol,visl)
  call psat(t, sat_pressure, ierr)
   call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
  !  print *,'visw  ',visl,tmp
  ! dif= 1.D-7 !effective diffusion coeff in liquid phase: check pcl
  

  
  diff(1:num_spec) = dif(1)
    !apply mixing rules
    ! should be careful with the consistance of mixing rules

!FEHM mixing
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol

! ideal mixing    
  !  den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

  den(1) = dw_mol
  h(1) = hw ! * xmol(1) + hg*xmol(2) 
  u(1) = h(1) - pw /den(1)* energyscale
  diff(1:num_spec) = dif(1)
  kvr(1) = kr(1)/visl
    !xlw = 1.D0
    !avgmw(1) = xlw*fmwh2o+(1.d0-xlw)*fmwco2
 ! no more gas phase
 !print *,'Trans-richards', x, t,p, satu
 ! note nphase=1, nspec = 1 + num_tracer, ndof = nspec +1  
 ! if assign nphase=2, it will still run, but gas phase properties will be  
 ! meaningless, wasting memory.
 
  nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
end subroutine pri_var_trans_Richards_ninc_2_2

 
subroutine pri_var_trans_Richards_ninc(x,iphase,energyscale,num_phase,num_spec, &
                                  ipckrtype,pckr_sir,pckr_lambda,pckr_alpha, &
                                  pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif, &
                                  var_node,itable,ierr,pref)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  implicit none
  
  integer :: num_phase,num_spec,num_pricomp
  integer :: size_var_use
  real*8 :: x(1:num_spec+1),energyscale
  real*8 :: var_node(1:2 + 7*num_phase + 2* num_phase*num_spec)
  real*8 :: dif(:)
  integer :: iphase, itable,ierr
  integer :: ipckrtype !, ithrmtype
     
    
  real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr 
  real*8 :: pref
  
 

  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
!  if ((num_phase == 1).and.( num_spec == 2)) then
    call pri_var_trans_Richards_ninc_2_2(x,iphase,energyscale,num_phase,num_spec, &
                                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha, &
                                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif, &
                                    var_node,itable,ierr, pref)
 !  else 
 !   print *, 'Wrong phase-specise combination. Stop.'
 !   stop
!  endif

end subroutine pri_var_trans_Richards_ninc   
  
  
subroutine pri_var_trans_Richards_winc(x,delx,iphase,energyscale,num_phase,num_spec,&
                                  ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                                  pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
                                 var_node,itable,ierr, pref)

  implicit none

  integer :: num_phase,num_spec
  integer :: size_var_use,size_var_node

  real*8 :: x(1:num_spec+1),delx(1:num_spec+1),energyscale, pref
  real*8 :: var_node(:)
  real*8 :: dif(:)
  integer ::iphase,itable,ierr
  integer :: ipckrtype !, ithrmtype
   
  integer :: n
  real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr 
  real*8 xx(1:num_spec+1)

  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  size_var_node = (num_spec+2)*size_var_use
  
  do n=1, num_spec+1
    xx = x
    xx(n) = x(n)+ delx(n)
  ! note: var_node here starts from 1 to grid%ndof*size_var_use
    call pri_var_trans_Richards_ninc(xx,iphase,energyscale,num_phase,num_spec,&
                                ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                                pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
                                var_node((n-1)*size_var_use+1:n*size_var_use), &
                                itable,ierr,pref)
  enddo

end subroutine pri_var_trans_Richards_winc
  
end module translator_Richards_module
