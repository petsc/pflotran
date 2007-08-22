 module translator_ths_module
 
  
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

    
 public  pri_var_trans_ths_ninc,pri_var_trans_ths_winc ,translator_ths_step_maxchange, &
         translator_ths_get_output,translator_check_phase_cond,translator_thsase_massbal
     
     
 real*8, private, parameter :: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                              fmwco2 = 44.0098D0, fmwnacl = 58.44277D0
 real*8, private, parameter :: eps=5D-7 , formeps=1D-2
 real*8, private, parameter ::yh2o_in_co2=1D-2   
 real*8, private, parameter :: rgasj   = 8.3143    ![J/K/mol]

 contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules


 subroutine translator_ths_massbal(grid)
 
  use pflow_gridtype_module
  
  implicit none
  
  type(pflowGrid) :: grid 
  
 
  integer, save :: icall
  integer :: ierr
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
  n2p=0
  nxc=0.; nyc=0.; nzc=0.D0; nzm=grid%z(grid%nmax); sm=0.; nsm=nzm; c0=0.D0

  do n = 1,grid%nlmax
    n0=(n-1)* grid%ndof
    index=(n-1)*size_var_node
    den=>var_p(index+3+grid%nphase: index+2+2*grid%nphase)
    sat=>var_p(index+2+1:index+2+grid%nphase)
    xmol=>var_p(index+2+7*grid%nphase+1:index+2+7*grid%nphase +&
         grid%nphase*grid%nspec)    
   
    pvol=volume_p(n)*porosity_p(n)         
    if(dabs(iphase_p(n)- 3.D0)<.25D0)then
      n2p=n2p+1
      x=grid%x(grid%nL2A(n)+1)
      y=grid%y(grid%nL2A(n)+1)
      z=grid%z(grid%nL2A(n)+1)
   
      if(z<nzm) nzm=z
      if(sm<sat(2))then
      !print *, n,grid%nL2A(n)+1,x,y,z,sat(2)
        sm=sat(2);nsm=z
      endif    
      c0=c0+sat(2)*pvol
      nxc=nxc + pvol*sat(2)*x
      nyc=nyc + pvol*sat(2)*y
      nzc=nzc + pvol*sat(2)*z
    endif

    do nc =1,grid%nspec
      do np=1,grid%nphase
        sum= sat(np)* xmol((np-1)*grid%nspec +nc)*den(np)
        tot(nc,np)= pvol*sum + tot(nc,np)
       !tot(0,np)=tot(0,np)+tot(nc,np)
       !tot(nc,0)=tot(nc,0)+tot(nc,np)
      enddo
    enddo
    nullify(sat, den,xmol) 
!   print *,nzc,c0
  enddo
 !  call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(grid%var,var_p,ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
 
 
  if(grid%commsize >1)then
    call MPI_REDUCE(n2p, n2p0,1,&
        MPI_INTEGER,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(c0, c00,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nxc, nxc0,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nyc, nyc0,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)

    call MPI_REDUCE(nzc, nzc0,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nzm, nzm0,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nsm, nsm0,1,&
        MPI_DOUBLE_PRECISION,MPI_MIN,0, PETSC_COMM_WORLD,ierr)
      
    do nc = 0,grid%nspec
      do np = 0,grid%nphase
        call MPI_REDUCE(tot(nc,np), tot0(nc,np),1,&
            MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
   
!       call MPI_BCAST(tot0,(grid%nphase+1)*(grid%nspec+1),&
!            MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
      enddo
    enddo
    if(grid%myrank==0) then
      tot = tot0; n2p=n2p0;nxc=nxc0;nyc=nyc0;nzc=nzc0;nzm=nzm0;nsm=nsm0;c0=c00 
    endif   
  endif 
  
  
  if(grid%myrank==0)then
    if(c0<1D-6) c0=1.D0
    nxc=nxc/c0; nyc=nyc/c0;nzc=nzc/c0
  
  
    write(*,'(" Total CO2: t= ",1pe12.4," dt= ",1pe12.4," liq:",1pe13.6,&
   &" gas:",1pe13.6, " tot:", 1p2e13.6, " [kmol]",1p3e13.6)') &
    grid%t/grid%tconv,grid%dt/grid%tconv,tot(2,1),tot(2,2),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,'(''# time   dt   totl   totg   tot n2p'')')
      icall = 1
    endif
!   write(13,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
! &  " gas:",1pe13.6," tot:",1p2e13.6," [kmol]")')&
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2)
    write(13,'(1p9e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,&
    tot(2,1),tot(2,2),tot(2,1)+tot(2,2),real(n2p),nzc,nzm,nsm 
  endif    
  
  
  
  
 end subroutine translator_ths_massbal



 subroutine translator_ths_get_output(grid)
  use pflow_gridtype_module
  type(pflowGrid), intent(inout) :: grid
  
  PetscScalar, pointer :: t_p(:),p_p(:),c_p(:),s_p(:),cc_p(:),var_P(:)
  integer n, index_var_begin ,jn, size_var_node
! PetscScalar, pointer :: p,t,satu(:),xmol(:)
  
  call VecGetArrayF90(grid%var, var_p, ierr)
  call VecGetArrayF90(grid%pressure, p_p, ierr)
  call VecGetArrayF90(grid%temp, t_p, ierr)
  call VecGetArrayF90(grid%xmol, c_p, ierr)
  call VecGetArrayF90(grid%sat, s_p, ierr)
  call VecGetArrayF90(grid%conc, cc_p, ierr)
  !print *,' translator_ths_get_output gotten pointers'
  
  size_var_node=(grid%ndof+1)*(2+7*grid%nphase +2*grid%nphase*grid%nspec)
  
  do n = 1, grid%nlmax
    index_var_begin= (n-1) * size_var_node
    jn = 1 + (n-1)*grid%nphase 
    
  p_p(jn) = var_p(index_var_begin + 2) - var_p(index_var_begin+5*grid%nphase+3)
  p_p(jn+1) = var_p(index_var_begin + 2) - var_p(index_var_begin+5*grid%nphase+4)
   
    t_p(n)=var_p(index_var_begin + 1)

    c_p(jn)=var_p(index_var_begin+7*grid%nphase+4)
  c_p(jn+1)= var_p(index_var_begin+7*grid%nphase+6)
    cc_p(n)=c_p(jn+1)
  
  s_p(jn)=var_p(index_var_begin + 3) 
    s_p(jn+1)=var_p(index_var_begin + 4) 
 enddo
 
  call VecRestoreArrayF90(grid%var, var_p, ierr)
  call VecRestoreArrayF90(grid%pressure, p_p, ierr)
  call VecRestoreArrayF90(grid%temp, t_p, ierr)
  call VecRestoreArrayF90(grid%xmol, c_p, ierr)
  call VecRestoreArrayF90(grid%sat, s_p, ierr)
  call VecRestoreArrayF90(grid%conc, cc_p, ierr)
 
! work only for 2 phases
 end subroutine translator_ths_get_output


 subroutine translator_ths_step_maxchange(grid)
   use pflow_gridtype_module
   type(pflowGrid), intent(inout) :: grid
  

  PetscScalar, pointer :: xx_p(:), yy_p(:), iphase_p(:),var_p(:),iphase_old_p(:)
  real*8 :: comp1,comp,cmp  
! real*8 :: dsm,dcm  
  real*8 :: dsm0,dcm0  
  integer n
! integer j

   call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
    call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)

  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)
  call VecGetArrayF90(grid%iphas_old, iphase_old_p,ierr)
  call VecGetArrayF90(grid%var, var_p, ierr)
  
  comp=0.D0;comp1=0.D0
  do  n=1, grid%nlmax
    n0=(n-1)*grid%ndof 
    if(int(iphase_p(n)) == int(iphase_old_p(n)))then
     cmp=dabs(xx_p(n0+3)-yy_p(n0+3))
     if(int(iphase_p(n))==1 .or.int(iphase_p(n))==2)then
       if(comp<cmp) comp=cmp
     
    endif   
       if(int(iphase_p(n))==3)then
       if(comp1<cmp) comp1=cmp
      
    endif   
    else
!  print *,'phase changed', n, iphase_p(n), iphase_old_p(n)

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
 end  subroutine translator_ths_step_maxchange


  

 subroutine pri_var_trans_ths_ninc_exec(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg, dif, var_node,itable,m_nacl,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    use water_eos_module
    use gas_eos_module  
    use pckr_module
    use co2eos_module
    use span_wagner_module


    implicit none
    integer :: num_phase,num_spec, itable, ierr
    integer :: size_var_use 
    real*8 x(1:num_spec+1),energyscale
    real*8, target:: var_node(:)
    integer ::iphase
    integer :: ipckrreg !, ithrmtype
    real*8 :: dif(:)

   
 !   integer size_var_node = (grid%ndof+1)*size_var_use

  real*8, pointer :: t ,p
  real*8, pointer :: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
  real*8, pointer :: xmol(:),satu(:),diff(:)
  integer ibase 
  
! real*8 p1,tmp,co2_phi,co2_poyn,stea,dstea_p,dstea_t,hstea_p,hstea_t,dstea
! real*8 pckr_swir,xla
  real*8 p2
  real*8 pw,dw_kg,dw_mol,hw,sat_pressure,visl,xphi,dco2
  real*8 dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  real*8 ug
  real*8 henry,m_nacl
  real*8 dsteamol,hstea
  real*8 kr(num_phase)
  real*8 err,vphi,xm_nacl,x_nacl
  real*8 temp
  
  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  !pckr_swir=pckr_sir(1)
  
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
  ibase=ibase+num_phase*num_spec;
  diff=>var_node(ibase:ibase+num_phase*num_spec-1)


   p = x(1)
   t = x(2)
   satu(1:num_phase-1)= x(3:grid%ndof)
  ! if(xmol(2)<0.D0) xmol(2)=0.D0
  ! if(xmol(2)>1.D0) xmol(2)=1.D0
  
  temp = 0D0
  do m = 1,num_phase-1
    temp = temp + satu(m)
  enddo
  satu(num_phase) =1.D0 -temp  
    
    
 
    call PSAT(t, sat_pressure, ierr)
 !   initial guess
   
    err=1.D0

!  print *, 'in 2 phase solver'
 
   !  p2=p*xmol(4)
  !call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)

  
  !do while(err>1E-8)
  
     p2=p!*xmol(4)
 !   p2=p*xmol(4)

 !   if(p2>=5d4)then
 !     !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
 !    
 !     call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
 !     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,itable)
 !     dg= dg / fmwco2
 !     fg= fg * 1.D6 
 !     hg= hg * fmwco2
 !   else      
      call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
 !     fg=p2
 !   endif

  !  xphi = fg/p2
!   xphi  = 1.d0 ! test-pcl
    
!    call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
  !  call Henry_duan_sun(p2*1D-5, t,  henry) ! henry = mol/kg/bars
   
  !   call Henry_duan_sun(t, p2 *1D-5, henry,xphi,m_nacl, m_nacl,sat_pressure*1D-5)
      
  !   henry= 1D0 / (fmwh2o *1D-3) / (henry*1D-5 )/xphi 
     
     !print *, "translator :; henry ", henry, xphi, p*.99/henry        
     !note: henry = H/phi
     
  !   tmp = (p-sat_pressure)/(henry-sat_pressure)
 !    err= dabs(tmp-xmol(4))
!   tmp=xmol(4)
 ! enddo
  
      
 xmol(:)=1.D0
 
  
   ! if(p2<5d4) call visco2(t,dg*fmwco2,visg)
!***************  Liquid phase properties **************************
!  avgmw(1)= xmol(1)* fmwh2o + xmol(2) * fmwco2 
!  avgmw(2)= xmol(3)* fmwh2o + xmol(4) * fmwco2 
  

  avgmw(1)= fmwh20
  avgmw(2)= fmwa


   ! pure water
  !  pw = p   
    if(num_phase>=2)then
      call pflow_pckr_noderiv(num_phase,ipckrreg,satu,pc,kr)
      pw=p !-pc(1)
     ! print *, num_phase,ipckrreg,satu,pc,kr
     end if

    
  !  call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
   !    call VISW(t,pw,sat_pressure,visl,tmp,tmp2,ierr)
   ! call VISW_FLO(t,dw_mol,visl)
   ! call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
   !  print *,'visw  ',visl,tmp
   ! dif= 1.D-7 !effective diffusion coeff in liquid phase: check pcl


    diff(:)= 0.D0
    
    !apply mixing rules
    ! should be careful with the consistance of mixing rules

!FEHM mixing
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
  !den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing
!   x_nacl = m_nacl/(m_nacl + 1D3/fmwh2o)
   
   ! xmol(1) = xh2o + xnacl
!   avgmw(1)= (xmol(1) - x_nacl) * fmwh2o + x_nacl * fmwnacl + xmol(2) * fmwco2
!   vphi=1D-6*(37.51D0 + t*(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))

  ! den(1)=dw_kg/(1D0-(fmwco2*1D-3-dw_kg*vphi)*xmol(2)/(avgmw(1)*1D-3))
  ! den(1)=den(1)/avgmw(1)
  
      
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 
 !   den(1)= 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1)=dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1)=den(1)/avgmw(1)
  pw=p
  call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
  call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
    den(1) = dw_mol
    h(1) = hw 
    u(1) = h(1) - pw /dw_mol* energyscale
    diff(1:num_spec) = 0.D0
    kvr(1) = kr(1)/visl
    !xlw = 1.D0
    !avgmw(1) = xlw*fmwh2o+(1.d0-xlw)*fmwco2

!*****************Gas phase properties 2*************************
 
    
  !  if(xmol(3)>eps)then
  !     call steameos(t,p,p2,dstea,dsteamol,dstea_p,dstea_t,&
  !          hstea,hstea_p,hstea_t,energyscale,ierr) 
!print *,'steameos', t,p,pa,p-pa,sat_pressure, dsteamol,hstea,dg,hg,dsteamol/xgw
  !     dsteamol=dsteamol
  !     hstea=hstea!/xgw
!       tmp=dsteamol
   ! else
 !       dsteamol=dg /xmol(4)*xmol(3)
 !       hstea=hg 
    !end if


!if((p*xgw)>sat_pressure)then 

 !call steameos(t,p,p-sat_pressure,dstea,dsteamol,dstea_p,dstea_t,&
 !                 hstea,hstea_p,hstea_t,energyscale,ierr) 
! dsteamol=tmp
 !   dstea_p=dg_p
 !   dstea_t=dg_t
!  hstea=tmp
!    hstea_p=hg_p
!    hstea_t=hg_t
!  Now just make the steam same as ideal gas   
 !endif
!**********************************************************************


  ! den(2)= dg*xmol(4) + dw_mol*xmol(3)
   den(2)= dg 

!   call visgas_noderiv(t,pa,p,den(2),visg)
!call visgas_noderiv(t,pa,p,den(2),visg)
! call visco2(t,dg ,visg)
!    visg=8.d-6
    kvr(2)=kr(2)/visg

  
!    den(2)=1.D0/( xga/dg + xgw/dsteamol)  !*c 
  !  den(2)=1.D0/( xga/dg + 1.D0/dsteamol)
       h(2)=  hg 
 !   h(2)= ( hg*xga  + hstea*xgw ) 
 !   h(2)= ( hg *xga + hstea ) 
    u(2)=  h(2)-p/den(2) * energyscale
    pc(2)=0
   !  print *,'gas phase nonder h::',t,p,h(2),u(2),hg,hstea
    diff(num_spec+1:num_spec*2)= 0.D0!kr(2)* dif(2) * 1.01325D5/p*((t+273.15)/273.15)**1.8D0 ! dirty part
     phase          1(water)          2(gas)            3(co2)
!   species      1    2    3       1    2    3        1   2    3
!   values                        1.0  1.0  1.0
!_________________________________________________________________

  
   nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
 end subroutine pri_var_trans_ths_ninc_2_2


 
 subroutine pri_var_trans_ths_ninc(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable,m_nacl,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  implicit none
  integer :: num_phase,num_spec
! integer :: num_pricomp
  integer :: size_var_use
  real*8 x(1:num_spec+1),energyscale
  real*8 var_node(1:2 + 7*num_phase + 2* num_phase*num_spec)
  real*8 :: dif(:), m_nacl
  integer ::iphase, itable,ierr
  integer :: ipckrreg !, ithrmtype
       
 ! real*8, optional :: phi_co2, den_co2  
  
  

  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
 ! if((num_phase == 2).and.( num_spec ==2)) then
    call pri_var_trans_ths_ninc_exec( x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable,m_nacl,ierr)
  !  else 
  !   print *, 'Wrong phase-specise combination. Stop.'
  !   stop
  ! endif
  ! print *, 'ninc: ',x, var_node
   
  end subroutine pri_var_trans_ths_ninc   
  
  
 subroutine pri_var_trans_ths_winc(x,delx,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable, m_nacl,ierr)
  integer :: num_phase,num_spec
  integer :: size_var_use,size_var_node
    

    real*8 x(1:num_spec+1),delx(1:num_spec+1),energyscale
    real*8 var_node(:),m_nacl
  real*8 :: dif(:)
  integer ::iphase,itable,ierr
  integer :: ipckrreg !, ithrmtype
   
   
    real*8 xx(1:num_spec+1)


    size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
    size_var_node = (num_spec+2)*size_var_use
  
   do n=1, num_spec+1
         xx=x;  xx(n)=x(n)+ delx(n)
  ! note: var_node here starts from 1 to grid%ndof*size_var_use
       call pri_var_trans_ths_ninc(xx,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg, dif,&
          var_node((n-1)*size_var_use+1:n*size_var_use),itable,m_nacl,ierr)
    enddo


 end subroutine pri_var_trans_ths_winc
 
end module   translator_ths_module
