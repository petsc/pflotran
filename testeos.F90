program main
use span_wagner_module

implicit double precision(a-h,o-z)

real*8 x(3), pckr_sir(2),dif(2), var_node(24)


   call initialize_span_wagner(0,0)
   p=2D7
   t=50.D0
   xmol=0.03
   
    x(1)=p;x(2)=t;x(3)=xmol
    dif(1)=1D-9;dif(2)=1D-6
    pckr_sir(1)=0.05;pckr_sir(2)=0.0


     
    call pri_var_trans_ninc_2_2(x,1,1D-6,2,2,&
                    4,pckr_sir,0.75,0.75,&
                    3.04,1D6,0.25,7.0, dif,&
					var_node,0,ierr,xphi)

print *, varnode

end




 subroutine pri_var_trans_ninc_2_2(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr, dif,&
					var_node,itable,ierr,xphi)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    use water_eos_module
    use gas_eos_module  
    use pckr_module
    use co2eos_module
    use span_wagner_module
  implicit none
real*8, parameter:: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                              fmwco2 = 44.0098D0
 real*8, parameter:: eps=5D-7
 real*8, parameter::yh2o_in_co2=1D-2   


  
    integer :: num_phase,num_spec, itable, ierr
	integer :: size_var_use 
    real*8 x(1:num_spec+1),energyscale
    real*8, target:: var_node(:)
	integer ::iphase
	integer :: ipckrtype !, ithrmtype
    real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
    real*8 :: dif(:)

   
 !   integer size_var_node = (grid%ndof+1)*size_var_use

    real*8, pointer :: t ,p
	real*8, pointer :: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
    real*8, pointer :: xmol(:),satu(:),diff(:)
    integer ibase 
  
    real*8 p1,p2,tmp
	real*8 pw,dw_kg,dw_mol,hw,sat_pressure,visl,xphi
	real*8 dg,dddt,dddp,fg, dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
	real*8 ug
	real*8 co2_phi, henry,co2_poyn
    real*8 stea,dsteamol,dstea_p,dstea_t, hstea,hstea_p,hstea_t,dstea
	real*8 kr(num_phase), pckr_swir
	real*8 err,xla,vphi

  
  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  pckr_swir=pckr_sir(1)
	
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

print *, size_var_use, iphase,x
	select case(iphase)
	

	case(1) ! only water phase
	 p = x(1)
	 t = x(2)
	 xmol(4)= x(3)
	! if(xmol(2)<0.D0) xmol(2)=0.D0
	! if(xmol(2)>1.D0) xmol(2)=1.D0
	if(xmol(4)<0.D0) print *,'tran:',iphase, x(1:3)
	if(xmol(4)>1.D0) print *,'tran:',iphase, x(1:3)
	 xmol(3)=1.D0 - xmol(4)
	! if(xmol(3)<0.D0)xmol(3)=0.D0
	 pc(1)=0.D0
	 pc(2)=0.D0
     satu(1)=1.D0
	 satu(2)= 0.D0
     !p2=p
	 ! if( p2>5d4)then
	 !  call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
     ! dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,itable)
     ! xphi=fg*1D6/p2
	 !else 
      !	  call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)
      !   call visco2(t,dg*fmwco2,visg)
	  !endif	 
	 
	 
   	! call Henry_CO2_noderiv(xla,tmp,t,p2, xphi, henry,co2_poyn)
	! print *,'Henry', xphi,p2, t, henry, co2_poyn  
		 
	  
     ! xmol(4)=henry*xmol(2)/p
	 ! xmol(3)=1.D0-xmol(4)
     ! if(xmol(1)<0.D0) xmol(1)=0.D0
	print *,'tran:'
	case(2) ! only supercritical CO2 Phase
		 
		 p = x(1)
    	 t = x(2)
	     xmol(4)= x(3)
		 if(xmol(4)<0.D0)print *,'tran:',iphase, x(1:3)
	     if(xmol(4)>1.D0) print *,'tran:',iphase, x(1:3)
		 
	     xmol(3)=1.D0 - xmol(4)
	!	 if(xmol(3)<0.D0)xmol(3)=0.D0
	 !    pc(1)=0.D0
		 pc(2)=0.D0
         satu(1)=eps
         satu(2)=1.D0
    
      
	
	case(3) ! water + gas 
   		 p = x(1)
    	 t = x(2)
	     satu(2)= x(3)
		 if(satu(2)<0.D0) print *,'tran:',iphase, x(1:3)
		 if(satu(2)>1.D0) print *,'tran:',iphase, x(1:3)
		 satu(1)=1.D0  - satu(2)
		 if( satu(2)<0.D0)satu(2)=0.D0
		 pc(2)=0.D0
     xmol(1)=1.D0; xmol(2)=0.D0; xmol(3)=yh2o_in_co2; xmol(4)=1.D0-xmol(3)

	
	end select	
		    
	  
    call PSAT(t, sat_pressure, ierr)

 !   initial guess
   
    err=1.D0

	print *, 'in 2 phase solver'
 
   !  p2=p*xmol(4)
  !call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)

  
  !do while(err>1E-8)  
     	  p2=p!*xmol(4)
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
     call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
  !   tmp = (p-sat_pressure)/(henry-sat_pressure)
 !    err= dabs(tmp-xmol(4))
!	 tmp=xmol(4)
 ! enddo
  
	    
 print *, 'out 2 phase solver'

    
 ! if(iphase /= 1)then
    xmol(3)=1.D0-xmol(4)
	if(xmol(3)<0.D0)xmol(3)=0.D0
	!if(xmol(3)<0.D0) xmol(3)=0.D0
	xmol(2)= p*xmol(4)/henry
    xmol(1)=1.D0-xmol(2)
	!if(xmol(1)<0.D0) xmol(1)=0.D0
  !endif
  if(p2<5d4) call visco2(t,dg*fmwco2,visg)
!***************  Liquid phase properties **************************
	avgmw(1)= xmol(1)* fmwh2o + xmol(2) * fmwco2 
    avgmw(2)= xmol(3)* fmwh2o + xmol(4) * fmwco2 
	

   ! pure water
    pw = p   
    if(num_phase>=2)then
     ! call pflow_pckr_noderiv(ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
     !                         pckr_m,pckr_pcmax,satu(2),pc,kr,pckr_betac,pckr_pwr)
      pw=p !-pc(1)

     end if

    
    call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
   !    call VISW(t,pw,sat_pressure,visl,tmp,tmp2,ierr)
   ! call VISW_FLO(t,dw_mol,visl)
    call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
  !  print *,'visw  ',visl,tmp
  ! dif= 1.D-7 !effective diffusion coeff in liquid phase: check pcl
	diff(1:num_spec)=dif(1)
    diff(num_spec +1: 2*num_spec)=dif(2) ! m^2/s @ 25C

    !apply mixing rules
    ! should be careful with the consistance of mixing rules
    
	!den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 
     vphi=1D-6*(37.51D0+t*(-9.585D-2+t*(8.74D-4-t*5.044D-7)))
	den(1)=dw_mol*fmwh2o/(1D0-(fmwco2*1D-3-dw_mol*fmwh2o*vphi)*xmol(2)/(avgmw(1)*1D-3))
	print *, den(1)
	den(1)=den(1)/avgmw(1)
	print *, den(1)
	h(1) = hw * xmol(1) + hg*xmol(2) 
    u(1) = h(1) - pw /den(1)* energyscale
    diff(1:num_spec) = dif(1)
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
        dsteamol=dg /xmol(4)*xmol(3)
        hstea=hg 
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
   den(2)= 1.D0/(xmol(4)/dg + xmol(3)/dw_mol)

!   call visgas_noderiv(t,pa,p,den(2),visg)
!call visgas_noderiv(t,pa,p,den(2),visg)
! call visco2(t,dg ,visg)
!    visg=8.d-6
    kvr(2)=kr(2)/visg

  
!    den(2)=1.D0/( xga/dg + xgw/dsteamol)  !*c 
  !  den(2)=1.D0/( xga/dg + 1.D0/dsteamol)
       h(2)=  hg *xmol(4) + hw*xmol(3) 
 !   h(2)= ( hg*xga  + hstea*xgw ) 
 !   h(2)= ( hg *xga + hstea ) 
    u(2)=  h(2)-p/den(2) * energyscale
    pc(2)=0
   !  print *,'gas phase nonder h::',t,p,h(2),u(2),hg,hstea
    diff(num_spec+1:num_spec*2)=dif(2) ! dirty part
    
 !if(t>=100.D0) print *, p,t,xga,xgw,dg,dsteamol,hg, hstea, h(2)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! remember: henry coeff is arranged in the vector as 
!    phase          1(water)          2(gas)            3(co2)
!   species      1    2    3       1    2    3        1   2    3
!   values                        1.0  1.0  1.0
!_________________________________________________________________

	
   nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
 end subroutine pri_var_trans_ninc_2_2

