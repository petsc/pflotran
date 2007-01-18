 module translator_ims_module
 
  
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
! phase                             Primary Variables			index
!	e								p, T, X(e,a), X(e,c)          1
!	g                               p, T, X(g,a), X(g,c)          2 
!   l								p, T, X(l,a), X(l,c)		  4	
!   eg								p, T,  S(g),  X(g,c)	      3 
!   el								p, T,  S(l),  X(l,c)          5
!   lg								p, T,  S(g),  X(g,c)          6
!   egl								p, T,  S(g),  S(l)            7
!**********************************************************************


! phase index 1.e; 2. l; 3. g
! within each phase component index : 1. H2O; 2. CO2; 3. Air

	  
 public  pri_var_trans_ims_ninc,pri_var_trans_ims_winc ,translator_ims_step_maxchange, &
         translator_ims_massbal		 
		 
		 
 real*8, private, parameter :: eps=5D-7 , formeps=5D-5, Rg = 8.3145D0

 contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules


 subroutine translator_ims_massbal(grid)
 use pflow_gridtype_module
  implicit none
  type(pflowGrid) :: grid 
  
 
  integer :: ierr,icall
  integer :: n,n0,nc,np,n2p,n2p0
  real*8 x,y,z,nzm,nzm0, nxc,nxc0,c0, c00,nyc,nyc0,nzc,nzc0,nsm,nsm0,sm 
  integer :: index, size_var_node
     
  PetscScalar, pointer ::  var_p(:),&
                           porosity_p(:), volume_p(:)
                           
   
  real*8 ::  pvol,sum
  real*8, pointer ::  den(:),sat(:)
 
  real*8 :: tot(0:grid%nphase), tot0(0:grid%nphase)  
  data icall/0/

  call VecGetArrayF90(grid%var,var_p,ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
 
  size_var_node=(grid%ndof+1)*(2 + 4*grid%nphase)
  tot=0.D0
  n2p=0
  nxc=0.; nyc=0.; nzc=0.D0; nzm=grid%z(grid%nmax); sm=0.; nsm=nzm; c0=0.D0

  do n = 1,grid%nlmax
    n0=(n-1)* grid%ndof
    index=(n-1)*size_var_node
    den=>var_p(index+3+grid%nphase: index+2+2*grid%nphase)
    sat=>var_p(index+2+1:index+2+grid%nphase)
    
   
    pvol=volume_p(n)*porosity_p(n)		     
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

	 
  
      do np=1,grid%nphase
        sum= sat(np)* den(np)
        tot(np)= pvol*sum + tot(np)
	    !tot(0,np)=tot(0,np)+tot(nc,np)
	    !tot(nc,0)=tot(nc,0)+tot(nc,np)
	  enddo
  
      nullify(sat, den) 
! print *,nzc,c0
  enddo
 !  call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  
  call VecRestoreArrayF90(grid%var,var_p,ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
 
 
 
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
	    

      do np = 0,grid%nphase
        call MPI_REDUCE(tot(np), tot0(np),1,&
	          MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
   
!       call MPI_BCAST(tot0,(grid%nphase+1)*(grid%nspec+1),&
!	          MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
      enddo

	if(grid%myrank==0) then
	   tot = tot0; n2p=n2p0;nxc=nxc0;nyc=nyc0;nzc=nzc0;nzm=nzm0;nsm=nsm0;c0=c00 
	endif   
  endif 
  
  
  if(grid%myrank==0)then
   if(c0<1D-6) c0=1.D0
	nxc=nxc/c0; nyc=nyc/c0;nzc=nzc/c0
	
	
    write(*,'(" Total CO2: liq:",1p, e13.6,&
   &" gas:",1p, e13.6, " tot:", 1p, 2e13.6, " [kg]",1p, 3e13.6)') &
   tot(1),tot(2),tot(1)+tot(2) !,nzc,nzm,nsm
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,*) '# time   dt   totl   totg   tot n2p'
      icall = 1
    endif
!   write(13,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
! &	" gas:",1pe13.6," tot:",1p2e13.6," [kmol]")')&
! & grid%t/grid%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2)
    write(13,'(1p9e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,&
 	tot(1),tot(2),tot(1)+tot(2),real(n2p),nzc,nzm,nsm 
  endif    
  
  
  
  
 end subroutine translator_ims_massbal


 subroutine translator_ims_step_maxchange(grid)
   use pflow_gridtype_module
   type(pflowGrid), intent(inout) :: grid
  

  PetscScalar, pointer :: xx_p(:), yy_p(:)
  real*8 :: dsm,dcm, comp1,comp, cmp  
  real*8 :: dsm0,dcm0  
  integer n, j

   call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
    
	do i=1, grid%nphase-1
	  call VecStrideNorm(grid%dxx,i,NORM_INFINITY,comp1,ierr)
      if(grid%dsmax< comp1) grid%dsmax=comp1     
    enddo

 end  subroutine translator_ims_step_maxchange
 
  
 
   
	   

 subroutine pri_var_trans_1(x,temp,energyscale,num_phase,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,&
					var_node,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    use water_eos_module
    use gas_eos_module  
    !use pckr_module
    !use co2eos_module
    !use span_wagner_module


    implicit none
    integer :: num_phase, ierr
	integer :: size_var_use 
    real*8 x(1:num_phase),energyscale,temp
    real*8, target:: var_node(:)
	integer ::iphase
	integer :: ipckrtype !, ithrmtype
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
      
 !   integer size_var_node = (grid%ndof+1)*size_var_use

    real*8, pointer :: t ,p
	real*8, pointer :: den(:),pc(:),kvr(:)
    real*8, pointer :: satu(:)
    integer ibase, i 
  
    real*8 p1,p2,tmp
	real*8 vis(num_phase)
	
	real*8 kr(num_phase), pckr_swir
	real*8 err, se

  
  size_var_use = 2 + 4*num_phase
  pckr_swir=pckr_sir(1)
	
	  
  ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;         p=>var_node(ibase)
  ibase=ibase+1;         satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
   

	 p = x(1)
	 t = temp ! kept the room for temperature for future density and viscosity functions
	  
 !    satu(1:num_phase-1)=x(2:num_phase)
	 tmp=0D0
	 do i=1, num_phase-1
	   
	   satu(i)=x(1+i)
	   if(satu(i)<0.D0) satu(i)=0.D0
	   if(satu(i)>1.D0) satu(i)=1.D0
	   tmp = tmp + satu(i)
	 enddo  
	     
	 satu(num_phase) = 1.D0 - tmp
     if(satu(num_phase)<0.D0) satu(num_phase)=0.D0
     if(satu(num_phase)>1.D0) satu(num_phase)=1.D0

	 
	    
	 pc(:)=0.D0 ! no capallary force considered now
	 kr(:)=1.D0
	 
     err=1.D0
     !call wateos_noderiv (t,p,p1,p2,tmp,energyscale,ierr) 
     !den(1)=p1
	 den(1)=1D3
	 vis(1)=5.494D-3 !, [Pa.S] water viscosity at 50C
	 if(num_phase>=2)then 
	 ! call ideal_gaseos_noderiv(p,t,energyscale, p1,p2, tmp)
	 ! den(2)=p1*44D0
      den(2)=30D0
	   vis(2)=14D-6    ! [Pa.S] CO2 vis
	 endif  
    !  if(num_phase>=3)then 
	 ! call ideal_gaseos_noderiv(p,t,energyscale, p1,p2, tmp)
	 ! den(2)=p1*44D0
     ! den(3)=.995D3
	 ! vis(3)=5.494D-2    ! [Pa.S] CO2 vis
	! endif  

     
     
     
    if(num_phase>=2)then
 !  Currently the relative permeability functions is exp functions of effective saturation.
        
        !tmp=0.D0; do i=1, num_phase; tmp = tmp + se(i);enddo
		
		 tmp =exp(1.D0) -1.D0
		do i=1, num_phase
		   kr(i)=0.D0
		   se =(satu(i)-pckr_sir(i))/(1.D0- pckr_sir(i))
		 !  if(se >=0.D0) kr(i) = (exp(se)-1.D0)/tmp 
		     kr(i)=se
		     kvr(i)=kr(i) / vis(i) 
		enddo
    end if

    
  	
   nullify(t, p, satu, den, pc,kvr)
 end subroutine pri_var_trans_1




 
 subroutine pri_var_trans_ims_ninc(x,temp,energyscale,num_phase,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,&
					var_node,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    implicit none
    integer :: num_phase
    integer :: size_var_use
	real*8 x(1:num_phase), temp, energyscale
    real*8 var_node(1:2 + 4*num_phase )
	integer :: ierr
	integer :: ipckrtype !, ithrmtype
     
	  
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr 
    !real*8, optional :: phi_co2, den_co2
	


    size_var_use = 2 + 4*num_phase 
    call pri_var_trans_1( x,temp, energyscale,num_phase,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,&
					var_node,ierr)
	
	end subroutine pri_var_trans_ims_ninc
!==================================================================================================



	
	
 subroutine pri_var_trans_IMS_winc(x,delx,temp,energyscale,num_phase,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,&
					var_node,ierr)
	integer :: num_phase
	integer :: size_var_use,size_var_node
    

    real*8 x(1:num_phase),delx(1:num_phase),energyscale, temp
    real*8 var_node(:)
	integer ::iphase,itable,ierr
	integer :: ipckrtype !, ithrmtype
   
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr 

    real*8 xx(1:num_phase)


    size_var_use = 2 + 4*num_phase 
	size_var_node = (num_phase)*size_var_use
	
	 do n=1, num_phase
         xx=x;  xx(n)=x(n)+ delx(n)
	! note: var_node here starts from 1 to grid%ndof*size_var_use
	    call pri_var_trans_ims_ninc(xx,temp,energyscale,num_phase,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,&
				var_node((n-1)*size_var_use+1:n*size_var_use),ierr)
	  enddo								
																										
					
 end subroutine pri_var_trans_ims_winc
	
end module 	translator_ims_module
