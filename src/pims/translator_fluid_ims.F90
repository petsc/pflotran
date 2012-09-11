 module translator_ims_module
 
  
 private 
	 
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

	  
 public  pri_var_trans_ims_ninc,pri_var_trans_ims_winc		 
		 
		 
 real*8, private, parameter :: eps=5D-7 , formeps=5D-5, Rg = 8.3145D0

 contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules

	   

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
    real*8 x(1:num_phase),energyscale,temp
    real*8, target:: var_node(:)
	integer ::iphase
	integer :: ipckrtype !, ithrmtype
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
      
 !   integergrid%size_var_node = (grid%ndof+1)*size_var_use

    real*8, pointer :: t ,p
	real*8, pointer :: den(:),pc(:),kvr(:)
    real*8, pointer :: satu(:)
    integer ibase, i 
  
    real*8 p1,p2,tmp
	real*8 vis(num_phase)
	
	real*8 kr(num_phase), pckr_swir
	real*8 err, se

  
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
     call wateos_noderiv (t,p,p1,p2,tmp,energyscale,ierr) 
     !den(1)=p1
	 den(1)=1D3
	 vis(1)=5.494D-3 !, [Pa.S] water viscosity at 50C
	 if(num_phase>=2)then 
	  call ideal_gaseos_noderiv(p,t,energyscale, p1,p2, tmp)
	  !den(2)=p1*44D0
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
