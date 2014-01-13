  module span_wagner_module

  ! module contains only 1 interface with other part, read as:

  ! co2_span_wagner(p,t,rho,dddt,dddp,fg,dfgdp,dfgdt,
  !                 eng,ent,dhdt,dhdp,visc,dvdt,dvdp)

  !  where units are assigned as:   
  !     P [MPa]       T [K]         
  !     rho [kg/m^3] Energy [MJ/kg] Enthalpy [MJ/kg] Vis [Pa s]

      implicit none
      save

      real*8, private :: n(42),ti(40),gamma(5),phic(8),c(40),d(40),a(8)
      real*8, private :: alpha(5),beta(8),delta(4),epsilon(5)
      real*8, private :: aco2(4),bco2(4),capa(5),capb(5),capc(5),capd(5)
      real*8, private :: av(0:4)

      real*8, private :: denc,tc,rg,pc
      real*8 ,private, allocatable:: co2_prop_spwag(:,:,:)
      real*8, private:: p,t,rhosav
	  public initialize_span_wagner, co2_span_wagner

      private

      contains

! ************************************************************************** !

subroutine initialize_span_wagner(itable,myrank)

      implicit none
      integer, optional:: itable
      
      real*8 pl,tl,tmp,tmp2
       
      integer iitable,i,j,myrank
      
      iitable=0
      if(present(itable)) iitable=itable

      denc = 467.6d0
      tc = 304.1282d0
      rg = 0.1889241d0
      pc = 7.3773d0

      av(0) = 0.235156d0
      av(1) = -.491266d0
      av(2) = 5.211155d-2
      av(3) = 5.347906d-2
      av(4) = -1.537102d-2

      a(1) = 8.37304456d0
      a(2) = -3.70454304d0
      a(3) = 2.5d0
      a(4) = 1.99427042d0
      a(5) = 0.62105248d0
      a(6) = 0.41195293d0
      a(7) = 1.04028922d0
      a(8) = 0.08327678d0

      phic(1) = 0.0d0
      phic(2) = 0.0d0
      phic(3) = 0.0d0
      phic(4) = 3.15163d0
      phic(5) = 6.1119d0
      phic(6) = 6.77708d0
      phic(7) = 11.32384d0
      phic(8) = 27.08792d0
      n(1) = 0.38856823203161d0
      n(2) = 0.2938547594274d1
      n(3) = -0.55867188534934d1
      n(4) = -0.76753199592477d0
      n(5) = 0.31729005580416d0
      n(6) = 0.54803315897767d0
      n(7) = 0.12279411220335d0
      n(8) = 0.2165896154322d1
      n(9) = 0.15841735109724d1
      n(10) = -0.23132705405503d0
      n(11) = 0.58116916431436d-1
      n(12) = -0.55369137205382d0
      n(13) = 0.48946615909422d0
      n(14) = -0.24275739843501d-1
      n(15) = 0.62494790501678d-1
      n(16) = -0.12175860225246d0
      n(17) = -0.37055685270086d0
      n(18) = -0.16775879700426d-1
      n(19) = -0.11960736637987d0
      n(20) = -0.45619362508778d-1
      n(21) = 0.35612789270346d-1
      n(22) = -0.74427727132052d-2
      n(23) = -0.17395704902432d-2
      n(24) = -0.21810121289527d-1
      n(25) = 0.24332166559236d-1
      n(26) = -0.37440133423463d-1
      n(27) = 0.14338715756878d0
      n(28) = -0.13491969083286d0
      n(29) = -0.2315122505348d-1
      n(30) = 0.12363125492901d-1
      n(31) = 0.2105832197294d-2
      n(32) = -0.33958519026368d-3
      n(33) = 0.55993651771592d-2
      n(34) = -0.30335118055646d-3
      n(35) = -0.2136548868832d3
      n(36) = 0.26641569149272d5
      n(37) = -0.24027212204557d5
      n(38) = -0.28341603423999d3
      n(39) = 0.21247284400179d3
      n(40) = -0.66642276540751d0
      n(41) = 0.72608632349897d0
      n(42) = 0.55068668612842d-1
      c(1) = 0
      c(2) = 0
      c(3) = 0
      c(4) = 0
      c(5) = 0
      c(6) = 0
      c(7) = 0
      c(8) = 1
      c(9) = 1
      c(10) = 1
      c(11) = 1
      c(12) = 1
      c(13) = 1
      c(14) = 1
      c(15) = 1
      c(16) = 1
      c(17) = 2
      c(18) = 2
      c(19) = 2
      c(20) = 2
      c(21) = 2
      c(22) = 2
      c(23) = 2
      c(24) = 3
      c(25) = 3
      c(26) = 3
      c(27) = 4
      c(28) = 4
      c(29) = 4
      c(30) = 4
      c(31) = 4
      c(32) = 4
      c(33) = 5
      c(34) = 6
      d(1) = 1
      d(2) = 1
      d(3) = 1
      d(4) = 1
      d(5) = 2
      d(6) = 2
      d(7) = 3
      d(8) = 1
      d(9) = 2
      d(10) = 4
      d(11) = 5
      d(12) = 5
      d(13) = 5
      d(14) = 6
      d(15) = 6
      d(16) = 6
      d(17) = 1
      d(18) = 1
      d(19) = 4
      d(20) = 4
      d(21) = 4
      d(22) = 7
      d(23) = 8
      d(24) = 2
      d(25) = 3
      d(26) = 3
      d(27) = 5
      d(28) = 5
      d(29) = 6
      d(30) = 7
      d(31) = 8
      d(32) = 10
      d(33) = 4
      d(34) = 8
      d(35) = 2
      d(36) = 2
      d(37) = 2
      d(38) = 3
      d(39) = 3
      ti(1) = 0.d0
      ti(2) = 0.75d0
      ti(3) = 1.0d0
      ti(4) = 2.0d0
      ti(5) = 0.75d0
      ti(6) = 2.0d0
      ti(7) = 0.75d0
      ti(8) = 1.5d0
      ti(9) = 1.5d0
      ti(10) = 2.5d0
      ti(11) = 0.0d0
      ti(12) = 1.5d0
      ti(13) = 2.0d0
      ti(14) = 0.0d0
      ti(15) = 1.0d0
      ti(16) = 2.0d0
      ti(17) = 3.0d0
      ti(18) = 6.0d0
      ti(19) = 3.0d0
      ti(20) = 6.0d0
      ti(21) = 8.0d0
      ti(22) = 6.0d0
      ti(23) = 0.0d0
      ti(24) = 7.0d0
      ti(25) = 12.0d0
      ti(26) = 16.0d0
      ti(27) = 22.0d0
      ti(28) = 24.0d0
      ti(29) = 16.0d0
      ti(30) = 24.0d0
      ti(31) = 8.0d0
      ti(32) = 2.0d0
      ti(33) = 28.0d0
      ti(34) = 14.0d0
      ti(35) = 1.0d0
      ti(36) = 0.0d0
      ti(37) = 1.0d0
      ti(38) = 3.0d0
      ti(39) = 3.0d0
      alpha(1) = 25d0
      alpha(2) = 25d0
      alpha(3) = 25d0
      alpha(4) = 15d0
      alpha(5) = 20d0
      beta(1) = 325d0
      beta(2) = 300d0
      beta(3) = 300d0
      beta(4) = 275d0
      beta(5) = 275d0
      beta(6) = 0.3d0
      beta(7) = 0.3d0
      beta(8) = 0.3d0
      gamma(1) = 1.16d0
      gamma(2) = 1.19d0
      gamma(3) = 1.19d0
      gamma(4) = 1.25d0
      gamma(5) = 1.22d0
      epsilon(1) = 1d0
      epsilon(2) = 1d0
      epsilon(3) = 1d0
      epsilon(4) = 1d0
      epsilon(5) = 1d0
      aco2(1)=3.5d0
      aco2(2)=3.5d0
      aco2(3)=3.0d0
      bco2(1)=0.875d0
      bco2(2)=0.925d0
      bco2(3)=0.875d0
      capa(1) = 0.7d0
      capa(2) = 0.7d0
      capa(3) = 0.7d0
      capb(1) = 0.3d0
      capb(2) = 0.3d0
      capb(3) = 1.0d0
      capc(1) = 10.d0
      capc(2) = 10.d0
      capc(3) = 12.5d0
      capd(1) = 275d0
      capd(2) = 275d0
      capd(3) = 275d0
	  
	  
  allocate( co2_prop_spwag(0:500,0:100,1:15))
  
  if(iitable == 1) then
      ! Generate the table, 3 index array
         !  1 p,  2  T
         !  3 rho 4 dddt, 5 dddp,
         !  6 fg,  7dfgdp,8 dfgdt
         !  9 eng
         !  10 ent, 11 dhdt,  12 dhdp,
         ! 13  visc, 14 dvdt, 15 dvdp
          
    if (myrank==0) print *,' Preparing Table...'
  ! allocate( co2_prop_spwag(0:1000,0:100,1:15))
    tmp2=0.D0    
       
    do i=0, 500
      tmp=tmp2
      do j=0,100
        pl= 0.01D0 + 0.5D0* real(i)
!       pl= 0.01D0 + 0.05D0* real(i)
        tl= 35.D0 + 273.15D0 + 2.5D0 * real(j)
        
        co2_prop_spwag(i,j, 1) = pl
        co2_prop_spwag(i,j, 2) = tl
        co2_prop_spwag(i,j, 3) = tmp 
        call co2_span_wagner(pl,tl,co2_prop_spwag(i,j,3),co2_prop_spwag(i,j,4),&
        co2_prop_spwag(i,j,5),co2_prop_spwag(i,j,6),co2_prop_spwag(i,j,7),&
        co2_prop_spwag(i,j,8),co2_prop_spwag(i,j,9),co2_prop_spwag(i,j,10),&
        co2_prop_spwag(i,j,11),co2_prop_spwag(i,j,12),co2_prop_spwag(i,j,13),&
        co2_prop_spwag(i,j,14),co2_prop_spwag(i,j,15))

        tmp=co2_prop_spwag(i,j, 3)
        if(j==0) tmp2= co2_prop_spwag(i,j, 3)
      ! print *, co2_prop_spwag(i,j,:)
      enddo
    enddo
  endif

end subroutine initialize_span_wagner      

! ************************************************************************** !

subroutine co2_span_wagner(pl,tl,rho,dddt,dddp,fg,dfgdp,dfgdt, &
      eng,ent,dhdt,dhdp,visc,dvdt,dvdp,itable)
      
  use co2_sw_rtsafe_module

  implicit none
      
      real*8 :: pl,tl,rho,eng,ent,dhdt,dhdp,dddt,dddp,visc,dvdt,dvdp
      real*8 :: fiot,fiott,fird,firt,firdt,firdd,firtt,ftau,fdel,del,tau,fir
!     real*8 :: vpartial, rho_h2o, Cs, xco2, rho_wco2
      real*8 :: fg,dfgdp,fg1,dfgdt
      real*8 :: rho1,rho2
      
!     integer :: it
      integer, optional :: itable
       
      integer iitable,i1,i2,j1,j2,i,isucc
      real*8 iindex,jindex,tmp,factor(1:4)
      
	  
      p=pl;t=tl;iitable=0
      if(present(itable)) iitable=itable

!     units are 
!     P : MPa
!     T : K
!     h(ent) : MJ/Kg
!     e(eng) : MJ/Kg
!     dhdt : MJ/kg/C
!     dhdp : MJ/Kg/MPa
!     rho : kg/m3
!     dddt : kg/m3/C
!     dddp : kg/m3/MPa
!     visc : Pa.s
!     dvdt: Pa.s/C
!     dvdp : Pa.s/Mpa

!     print *,'span_wag: ',p,t,tc,pc,rg

!     call co2_density (rho,it)
!     call co2_density0 (p,t,rho)

!*************************** Table Lookup *********************

  tablelook:  if(iitable==1)then

  isucc=1
! tmp = (p - 1.D-2) / 5.d-2; i1 = floor(tmp); i2 = ceiling(tmp); iindex=tmp 
  tmp = (p - 1.D-2) / 5d-1; i1 = floor(tmp); i2 = ceiling(tmp); iindex=tmp 
  tmp = (t - 308.15D0 ) / 2.5D0; j1 = floor(tmp); j2 = ceiling(tmp); jindex=tmp 

  if(iindex > 500.d0 .or. iindex < 0.d0 .or. jindex < 0.d0 .or. jindex > 100.d0) then
    print  *,' Out of Table Bounds: ', p,t,iindex,jindex
    isucc=0
  endif

  if(isucc>0)then
  	
	if(i1==i2 .and. j1==j2)then
	  factor(1)=1.D0
	  factor(2:4)=0.D0
    elseif(i2==i1) then 
	  factor(1)=-(jindex-j2)
	  factor(2)= 0.D0
	  factor(3)= (jindex-j1)
	  factor(4)= 0.D0
    elseif(j2==j1) then 
	  factor(1)=-(iindex-i2) 
	  factor(2)=(iindex-i1) 
	  factor(3)=0.0
	  factor(4)=0.D0
	else
	  factor(1)= (iindex-i2) * (jindex-j2) 
      factor(2)= -(iindex-i1) * (jindex-j2) 
      factor(3)= -(iindex-i2) * (jindex-j1) 
      factor(4)= (iindex-i1) * (jindex-j1) 
    endif

   
	  
		  
    i=1
    tmp = factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)
    if (dabs(tmp-p)>1D-10 ) then
      print *,' Error in intropolate::P',tmp,p,iindex,factor;isucc=0
    endif
   !print *, 'Table: P ',iindex,jindex, factor,i  
    i=i+1
    tmp = factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)
    if (dabs(tmp-t)>1D-10 ) then
      print *,' Error in intropolate:;T', tmp,t,jindex,factor; isucc=0
    endif
  endif 
  !print *, 'Table: T',iindex,jindex, factor,i 
  if(isucc==1)then       
    i=i+1
    rho = factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)
    i=i+1 
    dddt = factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dddp= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    fg= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dfgdp= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dfgdt= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    eng= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    ent= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dhdt= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dhdp= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    visc= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dvdt= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)

    i=i+1
    dvdp= factor(1)* co2_prop_spwag(i1,j1,i) + factor(2)* co2_prop_spwag(i2,j1,i) &
         + factor(3)* co2_prop_spwag(i1,j2,i) + factor(4)* co2_prop_spwag(i2,j2,i)
! print *,'sp:tab: ', p,t,rho,ent, factor
    return 
    endif     
  endif tablelook


!*************************** Solving EOS *********************
      call guess(rho1,rho2)
!     print *,'spanwag-guess: ',p,t,rho1,rho2

      rho = rtsafe(co2den,rho1,rho2,1.d-8)

      del = rho/denc
      tau = tc/t
      call phir(fir,del,tau)
      call dphiodtau(fiot,del,tau)
      call dphirddel(fird,del,tau)
      call dphirdtau(firt,del,tau)
      call dphirdddel(firdd,del,tau)
      call dphirdtautau(firtt,del,tau)
      call dphirddeldtau(firdt,del,tau)
      call dphiodtautau(fiott,del,tau)
      
      ent = rg*t*(1.0d0+tau*(fiot+firt)+del*fird) + 506.78d0
      eng = ent - p/rho*1.e3

!     dhdt = -rg*tau*tau*((-1.d0/(tau*tau))+fiott+firtt+ &
!         ((del*firdt)/tau)-((del*fird)/(tau*tau)))
		  
      dhdt = rg*((1.d0+del*(fird-tau*firdt))**2/(1.d0+del*(2.d0*fird+del*firdd))- &
	         tau*tau*(fiott+firtt))
	  
      dhdp = (1000.d0/denc)*((tau*firdt)+(del*firdd)+fird)/ &
          (1.d0+(del*del*firdd)+(2.d0*del*fird))

!     change units to MJ/kg from KJ/Kg
      ent = ent*1.0d-3
      eng = eng*1.0d-3
      dhdt = dhdt*1.0d-3
      dhdp = dhdp*1.d-3

      ftau = (del*del*firdt)-(1000.d0*p/(denc*rg*tc))
      fdel = (2.d0*del*fird)+(del*del*firdd)+1.d0
      dddt = (denc*tc/(t*t))*(ftau/fdel)
      dddp = 1000.d0/(rg*t*(1.d0+(del*del*firdd)+(2.d0*del*fird)))

!     fugacity
      fg1 = exp(fir+(del*fird)-log(1.d0+(del*fird)))
      fg = fg1*p
      dfgdp = (fird*(1.d0+(2.d0*del*fird)+(del*del*firdd))/ &
          (1.d0+(del*fird)))*dddp/denc
      dfgdp = fg1+(fg*dfgdp)
      
      dfgdt = (firt+(del*firt*fird)+(del*del*firdt*fird))/ &
          (1.d0+(del*fird))
      dfgdt = -fg*dfgdt*tc/(t*t)
      
!      write(*,*) "Molality?"
!      read(*,*) mol
!     call dissco2(p,t,mco2,fg1,mol)

!     density of aqueous CO2 solution
!     vpartial = 37.51d0-(9.585d-2*(t-273.15d0))+ &
!     (8.74d-4*((t-273.15d0)**2.d0))-(5.044d-7*((t-273.15d0)**3))

!      write(*,*) "Density of water?"
!      read(*,*) rho_h2o

!     Cs = mco2*rho_h2o
!     Xco2 = Cs/(Cs+(rho_h2o/18.d-3))
!     rho_wco2 = rho_h2o + (44.d-3*Cs) - (vpartial*1d-6*rho_h2o*Cs)

!      write(*,*) mco2, rho_wco2, xco2
!     Xco2 = .0292d0
!     rho_wco2 = rho_h2o + (1.96d2*Xco2) + (1.54d4*Xco2*Xco2)
!      write(*,*) rho_wco2

      call viscosity(pl,tl,rho,dddp,dddt,visc,dvdt,dvdp)

!     change units to Pa.s from microPa.s
      visc = visc*1d-6
      dvdt = dvdt*1d-6
      dvdp = dvdp*1d-6
      
end subroutine co2_span_wagner

! ************************************************************************** !

subroutine guess(lguess,uguess)

      implicit none

      real*8 :: ts,dum1,dum2,uguess,lguess,rhov,rhol
      real*8 :: a1,a2,a3,a4,a5,t1,t2,t3,t4,t5,tr

      if (p.le.pc .and. t.le.tc) then
         call vappr (ts,p,dum1,dum2,12)

!vapor density
            a1 = -1.7074879d0
            a2 = -0.8227467d0
            a3 = -4.6008549d0
            a4 = -10.111178d0
            a5 = -29.742252d0
            t1 = 0.34d0
            t2 = 0.5d0
            t3 = 1.d0
            t4 = 7.d0/3.d0
            t5 = 14.d0/3.d0
!           tr = 1.d0-(ts+273.15d0)/tc
            tr = 1.d0-ts/tc
            rhov =a1*(tr**t1)+a2*(tr**t2)+a3*(tr**t3)+ &
                    a4*(tr**t4)+a5*(tr**t5)
            rhov =denc*dexp(rhov)

!liquid density
            a1 = 1.9245108d0
            a2 = -0.62385555d0
            a3 = -0.32731127d0
            a4 = -0.39245142d0
            t1 = 0.34d0
            t2 = 0.5d0
            t3 = 5.d0/3.d0
            t4 = 11.d0/6.d0
            rhol =a1*(tr**t1)+a2*(tr**t2)+a3*(tr**t3)+ &
                    a4*(tr**t4)
            rhol =denc*dexp(rhol)

!           print *,'guess: ',p,t,ts,tc,denc,tr,rhov,rhol
 
           if (t.le.ts) then ! sub-critical liquid region
             if (p.lt.6.d0) then
               lguess = rhol+2.d0*(ts-t)
             else if (p.lt.7.d0) then
               lguess = rhol+4.d0*(ts-t)
             else
               lguess = rhol+6.d0*(ts-t)
             endif
             uguess = 2000.d0
           else             ! vapor region
              uguess = 1.d-2
              lguess = rhov
           endif
      else if (t.le.tc .and. p.gt.pc) then
        if (p.le.8.d0) then
          uguess = 2000.d0
          if (t.lt.275.d0) then
            lguess = 950.d0
          else if (t.le.290.d0) then
            lguess = 750.d0
          else
            lguess = 650.d0
          endif
        else if (p.le.9.5d0) then
          uguess = 2000.d0
          if (t.lt.275.d0) then
            lguess = 950.d0
          else if (t.le.290.d0) then
            lguess = 750.d0
          else
            lguess = 650.d0
          endif
        else if (p.le.100.d0) then
          uguess = 2000.d0
          if (t.lt.275.d0) then
            lguess = 950.d0
          else if (t.le.290.d0) then
            lguess = 750.d0
          else
            lguess = 650.d0
          endif
        else
          uguess = 2000.d0
          lguess = 850.d0
        endif
      else if (t.gt.tc .and. p.le.pc) then
         uguess = 2000.d0
         lguess = 1.d-2
      else if (p.gt.pc) then ! supercritical
         uguess = 2000.d0
         lguess = 50.d0
      endif

end subroutine guess

! ************************************************************************** !

subroutine co2den(den,f,df)
     
	  IMPLICIT NONE
      real*8 :: den,tau1,del1
      real*8 :: f1,df1,f,df

      tau1 = tc/t
      del1 = den/denc

      call dphirddel(f1,del1,tau1)
      call dphirdddel(df1,del1,tau1)
      f = del1*(1.d0+del1*f1)-(p/(denc*0.001d0*rg*t))
      df = 1.d0+(2.d0*del1*f1)+(del1*del1*df1)

      return
end subroutine co2den

! ************************************************************************** !

double precision function psi(i,del2,tau2)
      implicit none
!     real*8 :: psi
!     real*8 aco2(4),bco2(4),capa(5),capb(5),capc(5),capd(5)
      real*8 del2,tau2
      integer i
!     common/params3/aco2,bco2,capa,capb,capc,capd

      psi=-(capc(i)*((del2-1.d0)**2.d0))-(capd(i)*((tau2-1.d0)**2.d0))
      psi=exp(psi)

      end function psi

! ************************************************************************** !

      subroutine dphiodtau(dr,del2,tau2)
      implicit none
      integer :: i
      real*8 :: del2,tau2,dr
      real*8 :: ideal_helm, derti_helm

      ideal_helm = log(del2)+a(1)+a(2)*tau2+a(3)*log(tau2)
      do i = 4, 8
        ideal_helm = ideal_helm+(a(i)*log(1.d0-exp(-phic(i)*tau2)))
      enddo

!     Table 34 derivative of ideal helmholtz wrt tau

      derti_helm = a(2) + (a(3)/tau2)
      do i = 4, 8
       derti_helm = derti_helm+a(i)*phic(i)*(1.d0/(1.d0-exp(-phic(i)*tau2))-1.d0)
      enddo
      dr = derti_helm

end subroutine dphiodtau

! ************************************************************************** !

subroutine dphiodtautau(dr,del2,tau2)

      implicit none
      
!     Span & Wagner (1996) Table 32

      integer :: i
      real*8 :: del2,tau2,dr,dihelm_dtautau

!     Table 34 double derivative of ideal helmholtz wrt tau

      dihelm_dtautau = -a(3)/(tau2*tau2)
      do i = 4, 8
        dihelm_dtautau = dihelm_dtautau-a(i)*(phic(i)**2)* &
        (1.d0-exp(-phic(i)*tau2))**(-2)*exp(-phic(i)*tau2)
      enddo
      dr = dihelm_dtautau

end subroutine dphiodtautau

! ************************************************************************** !

subroutine phir(r_helm,del2,tau2)

!     residual helmholtz energy: Span & Wagner (1996), p. 1544, eq. (6.5)
!     Table 32

      implicit none
      
      real*8 :: del2,tau2,r_helm,psi1
      real*8 :: x1,x2,x3,x4,x5,x6,x7,x8,x10
      real*8 :: e1,e2,e3,e4,e5,e6
      real*8 :: tsqr,t1,t2,t3,t6,t7,t8,t12,t14,t16,t22,t24,t28
      real*8 :: capdel1
      integer :: i

!     equation 5.3 for residual part of Helmholtz function

      r_helm = 0.0d0

!     do i = 1, 7
!       r_helm = r_helm+n(i)*(del2**d(i))*(tau2**ti(i))
!     enddo

      x1 = del2
      x2 = del2*del2
      x3 = del2*x2
      x4 = x2*x2
      x5 = del2*x4
      x6 = x2*x4
      x7 = del2*x6
      x8 = x2*x6
      x10 = x2*x8

 !    r_helm = n(1)*del2 + &
 !             n(2)*del2*tau2**0.75d0 + &
 !             n(3)*del2*tau2 + &
 !             n(4)*del2*tau2**2 + &
 !             n(5)*x2*tau2**0.75d0 + &
 !             n(6)*x2*tau2**2 + &
 !             n(7)*x3*tau2**0.75d0

      r_helm = (n(1) + &
               n(2)*tau2**0.75d0 + &
               n(3)*tau2 + &
               n(4)*tau2**2 + &
              (n(5)*tau2**0.75d0 + &
               n(6)*tau2**2 + &
               n(7)*del2*tau2**0.75d0)*del2)*del2

!     do i = 8, 34
!       r_helm = r_helm+n(i)*del2**d(i)*tau2**ti(i)*exp(-del2**c(i))
!     enddo

      e1 = exp(-del2)
      e2 = exp(-del2**2)
      e3 = exp(-del2**3)
      e4 = exp(-del2**4)
      e5 = exp(-del2**5)
      e6 = exp(-del2**6)

      tsqr = sqrt(tau2)
      t1 = tau2
      t2 = tau2*tau2
      t3 = tau2*t2
      t6 = t3*t3
      t7 = t1*t6
      t8 = t6*t2
      t12= t6*t6
      t14= t12*t2
      t16= t14*t2
      t22= t16*t6
      t24= t12*t12
      t28= t22*t6

      r_helm = r_helm + (((n(8)*del2 + &
                      n(9)*x2)*t1 + &
                      n(10)*x4*t2)*tsqr + &
                      n(11)*x5 + &
                      n(12)*x5*t1*tsqr + &
                      n(13)*x5*t2 + &
                      n(14)*x6 + &
                      n(15)*x6*tau2 + &
                      n(16)*x6*t2)*e1 + &
                     (n(17)*del2*t3 + &
                      n(18)*del2*t6 + &
                      n(19)*x4*t3 + &
                      n(20)*x4*t6 + &
                      n(21)*x4*t8 + &
                      n(22)*x7*t6 + &
                      n(23)*x8)*e2 + &
                     (n(24)*x2*t7 + &
                      n(25)*x3*t12 + &
                      n(26)*x3*t16)*e3 + &
                     (n(27)*x5*t22 + &
                      n(28)*x5*t24 + &
                      n(29)*x6*t16 + &
                      n(30)*x7*t24 + &
                      n(31)*x8*t8 + &
                      n(32)*x10*t2)*e4 + &
                      n(33)*x4*t28*e5 + &
                      n(34)*x8*t14*e6 

!     do i = 35, 39
!       r_helm = r_helm+n(i)*(del2**d(i))*(tau2**ti(i))* &
!       exp((-alpha(i-34)*((del2-epsilon(i-34))**2))-(beta(i-34)* &
!       ((tau2-gamma(i-34))**2)))
!     enddo

        r_helm = r_helm + (n(35)*tau2* &
        exp(-alpha(1)*(del2-epsilon(1))**2-beta(1)* &
        (tau2-gamma(1))**2) + &

        n(36)* &
        exp(-alpha(2)*(del2-epsilon(2))**2-beta(2)* &
        (tau2-gamma(2))**2) + &

        n(37)*tau2* &
        exp(-alpha(3)*(del2-epsilon(3))**2-beta(3)* &
        (tau2-gamma(3))**2) + &

        (n(38)* &
        exp(-alpha(4)*(del2-epsilon(4))**2-beta(4)* &
        (tau2-gamma(4))**2) + &

        n(39)* &
        exp(-alpha(5)*(del2-epsilon(5))**2-beta(5)* &
        (tau2-gamma(5))**2))*x1*t3)*x2

      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        r_helm=r_helm+n(i)*(capdel1**bco2(i-39))*del2*psi1
      enddo

      return
end subroutine phir

! ************************************************************************** !

subroutine dphirddel(dr,del2,tau2)

!     Span & Wagner (1996) Table 32

      implicit none
      integer :: i
      real*8 :: del2,tau2,dr,derdr_helm
      real*8 :: psi1,capdel1,dsidd,ddelbdd
      real*8 :: d1,d2,d3,d4,d5,d6,d7,d9
      real*8 :: e1,e2,e3,e4,e5,e6
      real*8 :: tsqr,t1,t2,t3,t6,t7,t8,t12,t14,t16,t22,t24,t28

!     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
!     do i = 1, 7
!       derdr_helm=derdr_helm+n(i)*d(i)*del2**(d(i)-1.d0)*tau2**ti(i)
!     enddo

      tsqr = sqrt(tau2)
      t1 = tau2
      t2 = tau2*tau2
      t3 = tau2*t2
      t6 = t3*t3
      t7 = t1*t6
      t8 = t6*t2
      t12= t6*t6
      t14= t12*t2
      t16= t14*t2
      t22= t16*t6
      t24= t12*t12
      t28= t22*t6

        derdr_helm=derdr_helm &
        + n(1) &
        + n(2)*tau2**0.75d0 &
        + n(3)*tau2 &
        + n(4)*t2 &
        + (n(5)*2.d0*tau2**0.75d0 &
        + n(6)*2.d0*t2 &
        + n(7)*3.d0*del2*tau2**0.75d0)*del2

        d1 = del2
        d2 = d1*d1
        d3 = d1*d2
        d4 = d2*d2
        d5 = d3*d2
        d6 = d3*d3
        d7 = d1*d6
        d9 = d2*d7

      e1 = exp(-del2)
      e2 = exp(-del2**2)
      e3 = exp(-del2**3)
      e4 = exp(-del2**4)
      e5 = exp(-del2**5)
      e6 = exp(-del2**6)

!     do i = 8, 34
!       derdr_helm=derdr_helm + n(i) * del2**(d(i)-1.d0) * tau2**ti(i) * &
!       exp(-del2**c(i)) * (d(i) - c(i) * del2**c(i))
!     enddo

        derdr_helm=derdr_helm + &
        (n(8)*t1*tsqr*(1.d0-del2) + &
        n(9)*del2*t1*tsqr*(2.d0-del2) + &
        n(10)*d3*t2*tsqr*(4.d0-del2) + &
        n(11)*d4*(5.d0-del2) + &
        n(12)*d4*t1*tsqr*(5.d0-del2) + &
        n(13)*d4*t2*(5.d0-del2) + &
        n(14)*d5*(6.d0-del2) + &
        n(15)*d5*tau2*(6.d0-del2) + &
        n(16)*d5*t2*(6.d0-del2))*e1 + &
        
        (n(17)*t3*(1.d0-2.d0*d2) + &
        n(18)*t6*(1.d0-2.d0*d2) + &
        n(19)*d3*t3*(4.d0-2.d0*d2) + &
        n(20)*d3*t6*(4.d0-2.d0*d2) + &
        n(21)*d3*t8*(4.d0-2.d0*d2) + &
        n(22)*d6*t6*(7.d0-2.d0*d2) + &
        n(23)*d7*(8.d0-2.d0*d2))*e2 + &
        
        (n(24)*del2*t7*(2.d0-3.d0*d3) + &
        n(25)*d2*t12*(3.d0-3.d0*d3) + &
        n(26)*d2*t16*(3.d0-3.d0*d3))*e3 + &
        
        (n(27)*d4*t22*(5.d0-4.d0*d4) + &
        n(28)*d4*t24*(5.d0-4.d0*d4) + &
        n(29)*d5*t16*(6.d0-4.d0*d4) + &
        n(30)*d6*t24*(7.d0-4.d0*d4) + &
        n(31)*d7*t8*(8.d0-4.d0*d4) + &
        n(32)*d9*t2*(10.d0-4.d0*d4))*e4 + &
        
        n(33)*d3*t28*e5*(4.d0-5.d0*d5) + &
        n(34)*d7*t14*e6*(8.d0-6.d0*d6)

!     do i = 35, 39
!       derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))* &
!       exp((-alpha(i-34)*((del2-epsilon(i-34))**2.0d0))-(beta(i-34)* &
!       ((tau2-gamma(i-34))**2.0d0)))* &
!       ((d(i)/del2)-(2.0d0*alpha(i-34)*(del2-epsilon(i-34)))))
!     enddo

!     do i = 35, 39
!       derdr_helm = derdr_helm + n(i)*del2**d(i)*tau2**ti(i)* &
!       exp(-alpha(i-34)*(del2-epsilon(i-34))**2 - beta(i-34)* &
!       (tau2-gamma(i-34))**2)*(d(i)/del2-2.d0*alpha(i-34)*(del2-epsilon(i-34)))
!     enddo

      derdr_helm = derdr_helm &
      + n(35)*d2*tau2* &
        exp(-alpha(1)*(del2-epsilon(1))**2 - beta(1)* &
        (tau2-gamma(1))**2)*(2.d0/del2-2.d0*alpha(1)*(del2-epsilon(1))) &

      + n(36)*d2* &
        exp(-alpha(2)*(del2-epsilon(2))**2 - beta(2)* &
        (tau2-gamma(2))**2)*(2.d0/del2-2.d0*alpha(2)*(del2-epsilon(2))) &

      + n(37)*d2*tau2* &
        exp(-alpha(3)*(del2-epsilon(3))**2 - beta(3)* &
        (tau2-gamma(3))**2)*(2.d0/del2-2.d0*alpha(3)*(del2-epsilon(3))) &

      + n(38)*d3*t3* &
        exp(-alpha(4)*(del2-epsilon(4))**2 - beta(4)* &
        (tau2-gamma(4))**2)*(3.d0/del2-2.d0*alpha(4)*(del2-epsilon(4))) &

      + n(39)*d3*t3* &
        exp(-alpha(5)*(del2-epsilon(5))**2 - beta(5)* &
        (tau2-gamma(5))**2)*(3.d0/del2-2.d0*alpha(5)*(del2-epsilon(5)))
      

      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        dsidd=dpsiddel(i-39,del2,tau2)
        ddelbdd=ddelbiddel(i-39,del2,tau2)
        derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i-39))*(psi1+ &
        (del2*dsidd)))+(del2*psi1*ddelbdd)))
      enddo

      dr = derdr_helm
      
end subroutine dphirddel

! ************************************************************************** !

subroutine dphirdddel(dpdd,del2,tau2)
      implicit none
      integer :: i
      real*8 :: del2,tau2,derdr_helm
      real*8 :: dpdd,psi1,capdel1
      real*8 :: dsidd,d2sidd,ddelbdd
      real*8 :: d2delbdd
      real*8 :: d1,d2,d3,d4,d5,d6,d7,d8,d9
      real*8 :: e1,e2,e3,e4,e5,e6
      real*8 :: tsqr,t1,t2,t3,t6,t7,t8,t12,t14,t16,t22,t24,t28

!     table 36, derivative of residual helmholtz wrt delta

		d1 = del2
		d2 = d1*d1
		d3 = d1*d2
		d4 = d2*d2
		d5 = d3*d2
		d6 = d3*d3
		d7 = d1*d6
		d8 = d2*d6
		d9 = d2*d7

      e1 = exp(-del2)
      e2 = exp(-del2**2)
      e3 = exp(-del2**3)
      e4 = exp(-del2**4)
      e5 = exp(-del2**5)
      e6 = exp(-del2**6)

      tsqr = sqrt(tau2)
      t1 = tau2
      t2 = tau2*tau2
      t3 = tau2*t2
      t6 = t3*t3
      t7 = t1*t6
      t8 = t6*t2
      t12= t6*t6
      t14= t12*t2
      t16= t14*t2
      t22= t16*t6
      t24= t12*t12
      t28= t22*t6

      derdr_helm = 0.0d0
!     do i = 1, 7
!       derdr_helm=derdr_helm+n(i)*d(i)*(d(i)-1.d0)* &
!       del2**(d(i)-2.d0)*tau2**ti(i)
!     enddo

        derdr_helm=derdr_helm &
!       + n(1)*d(i)*(d(i)-1.d0)*del2**(d(i)-2.d0)*tau2**ti(i) &
!       + n(2)*d(i)*(d(i)-1.d0)*del2**(d(i)-2.d0)*tau2**ti(i) &
!       + n(3)*d(i)*(d(i)-1.d0)*del2**(d(i)-2.d0)*tau2**ti(i) &
!       + n(4)*d(i)*(d(i)-1.d0)*del2**(d(i)-2.d0)*tau2**ti(i) &
        + n(5)*2.d0*tau2**ti(5) &
        + n(6)*2.d0*tau2**2 &
        + n(7)*6.d0*del2*tau2**ti(7)

!     do i = 8, 34
!        derdr_helm=derdr_helm+(n(i)*(del2**(d(i)-2))*(tau2**ti(i))*exp &
!        (-del2**c(i))*(((d(i)-(c(i)*(del2**c(i))))* &
!        (d(i)-1.d0-(c(i)*(del2**c(i)))))-((c(i)**2.d0)* &
!        (del2**c(i)))))
!     enddo
      
!     do i = 8, 34
!       derdr_helm=derdr_helm+n(i)*del2**(d(i)-2)*tau2**ti(i)*exp(-del2**c(i))* &
!       (((d(i)-c(i)*del2**c(i))* &
!       (d(i)-1.d0-c(i)*del2**c(i)))-(c(i)**2*del2**c(i)))
!     enddo

        derdr_helm=derdr_helm+ &
        n(8)*del2**(-1)*t1*tsqr*e1* &
        (((1.d0-del2)*(-del2))-del2) + &

        n(9)*t1*tsqr*e1* &
        (((2.d0-del2)*(1.d0-del2))-del2) + &

        n(10)*d2*t2*tsqr*e1* &
        (((4.d0-del2)*(3.d0-del2))-del2) + &

        n(11)*d3*e1* &
        (((5.d0-del2)*(4.d0-del2))-del2) + &

        n(12)*d3*t1*tsqr*e1* &
        (((5.d0-del2)*(4.d0-del2))-del2) + &

        n(13)*d3*t2*e1* &
        (((5.d0-del2)*(4.d0-del2))-del2) + &

        n(14)*d4*e1* &
        (((6.d0-del2)*(5.d0-del2))-del2) + &

        n(15)*d4*tau2*e1* &
        (((6.d0-del2)*(5.d0-del2))-del2) + &

        n(16)*d4*t2*e1* &
        (((6.d0-del2)*(5.d0-del2))-del2) + &

        n(17)*del2**(-1)*t3*exp(-del2**2)* &
        (((1.d0-2.d0*del2**2)* &
        (-2.d0*del2**2))-(4.d0*del2**2)) + &

        n(18)*del2**(-1)*t6*e2* &
        (((1.d0-2.d0*d2)* &
        (-2.d0*d2))-(4.d0*d2)) + &

        n(19)*d2*t3*e2* &
        (((4.d0-2.d0*d2)* &
        (3.d0-2.d0*d2))-(4.d0*d2)) + &

        n(20)*d2*t6*e2* &
        (((4.d0-2.d0*d2)* &
        (3.d0-2.d0*d2))-(4.d0*d2)) + &

        n(21)*d2*t8*e2* &
        (((4.d0-2.d0*d2)* &
        (3.d0-2.d0*d2))-(4.d0*d2)) + &

        n(22)*d5*t6*e2* &
        (((7.d0-2.d0*d2)* &
        (6.d0-2.d0*d2))-(4.d0*d2)) + &

        n(23)*d6*e2* &
        (((8.d0-2.d0*d2)* &
        (7.d0-2.d0*d2))-(4.d0*d2)) + &

        n(24)*t7*e3* &  ! ***check***
        (((2.d0-3.d0*d3)* &
        (2.d0-3.d0*d3))-(9.d0*d3)) + &

        n(25)*del2*t12*e3* &
        (((2.d0-3.d0*d3)* &
        (1.d0-3.d0*d3))-(9.d0*d3)) + &

        n(26)*del2*t16*e3* &
        (((3.d0-3.d0*d3)* &
        (2.d0-3.d0*d3))-(9.d0*d3)) + &

        n(27)*d3*t22*e4* &
        (((5.d0-4.d0*d4)* &
        (4.d0-4.d0*d4))-(16.d0*d4)) + &

        n(28)*d3*t24*e4* &
        (((5.d0-4.d0*d4)* &
        (4.d0-4.d0*d4))-(16.d0*d4)) + &

        n(29)*d4*t16*e4* &
        (((6.d0-4.d0*d4)* &
        (5.d0-4.d0*del2**4))-(16.d0*del2**4)) + &

        n(30)*del2**5*t24*exp(-del2**4)* &
        (((7.d0-4.d0*del2**4)* &
        (6.d0-4.d0*del2**4))-(16.d0*del2**4)) + &

        n(31)*d6*t8*e4* &
        (((8.d0-4.d0*d4)* &
        (7.d0-4.d0*d4))-(16.d0*d4)) + &

        n(32)*d8*t2*e4* &
        (((10.d0-4.d0*d4)* &
        (9.d0-4.d0*d4))-(16.d0*d4)) + &

        n(33)*d2*t28*e5* &
        (((4.d0-5.d0*d5)* &
        (3.d0-5.d0*d5))-(25.d0*d5)) + &

        n(34)*d6*t14*e6* &
        (((8.d0-6.d0*d6)* &
        (7.d0-6.d0*d6))-(36.d0*d6))

      
      do i = 35, 39
        derdr_helm = derdr_helm+(n(i)*(tau2**ti(i))* &
        exp((-alpha(i-34)*((del2-epsilon(i-34))**2))-(beta(i-34)* &
        ((tau2-gamma(i-34))**2)))* &
        ((-2.d0*alpha(i-34)*(del2**d(i)))+(4.d0*(alpha(i-34)**2)* &
        (del2**d(i))*((del2-epsilon(i-34))**2))-(4.d0*d(i)* &
        alpha(i-34)*(del2**(d(i)-1))*(del2-epsilon(i-34)))+ &
        (d(i)*(d(i)-1.d0)*(del2**(d(i)-2.d0)))))
      enddo
      
      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        dsidd=dpsiddel(i-39,del2,tau2)
        d2sidd=d2psiddel2(i-39,del2,tau2)
        ddelbdd=ddelbiddel(i-39,del2,tau2)
        d2delbdd=d2delbiddel2(i-39,del2,tau2)
        derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i-39))*((2.d0* &
        dsidd)+(del2*d2sidd)))+(2.d0*ddelbdd*(psi1+(del2*dsidd))) &
        +(d2delbdd*del2*psi1)))
      enddo
      
      dpdd = derdr_helm
      
end subroutine dphirdddel

! ************************************************************************** !

subroutine dphirdtau(dpdtau,del2,tau2)

      implicit none
      integer :: i
      real*8 :: del2,tau2,derdr_helm
      real*8 :: dpdtau
      real*8 :: psi1,dbidt,capdel1,dsidt

!     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)* &
        (tau2**(ti(i)-1))*(del2**d(i)))
      enddo
      
      do i = 8, 34
        derdr_helm=derdr_helm+(n(i)*ti(i)*(del2**d(i))* &
        (tau2**(ti(i)-1))*exp(-del2**c(i)))
      enddo
      
      do i = 35, 39
        derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))* &
        exp((-alpha(i-34)*((del2-epsilon(i-34))**2.0d0))-(beta(i-34)* &
        ((tau2-gamma(i-34))**2.0d0)))* &
        ((ti(i)/tau2)-(2.0d0*beta(i-34)*(tau2-gamma(i-34)))))
      enddo
      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        dbidt=ddelbidtau(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        dsidt=dpsidtau(i-39,del2,tau2)
        derdr_helm=derdr_helm+(n(i)*del2*((dbidt*psi1)+((capdel1** &
        bco2(i-39))*dsidt)))
      enddo
      dpdtau = derdr_helm
      
end subroutine dphirdtau

! ************************************************************************** !

subroutine dphirdtautau(dpdtt,del2,tau2)

      implicit none
      integer :: i
      real*8 :: del2,tau2,derdr_helm
      real*8 :: dpdtt
      real*8 :: psi1
      real*8 d2bidtt,capdel1,dbidt,dsidt,d2sidtt

!     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)*(ti(i)-1.d0)* &
        (tau2**(ti(i)-2.d0))*(del2**d(i)))
      enddo
      
      do i = 8, 34
        derdr_helm=derdr_helm+(n(i)*ti(i)*(ti(i)-1.d0)*(del2**d(i))* &
        (tau2**(ti(i)-2))*exp(-del2**c(i)))
      enddo
      
      do i = 35, 39
        derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))* &
        exp((-alpha(i-34)*((del2-epsilon(i-34))**2.0d0))-(beta(i-34)* &
        ((tau2-gamma(i-34))**2.0d0)))*( &
        (((ti(i)/tau2)-(2.0d0*beta(i-34)*(tau2-gamma(i-34))))**2.d0) &
        -(ti(i)/(tau2*tau2))-(2.d0*beta(i-34))))
      enddo
      
      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        d2bidtt=d2delbidtau2(i-39,del2,tau2)
        dbidt=ddelbidtau(i-39,del2,tau2)
        dsidt=dpsidtau(i-39,del2,tau2)
        d2sidtt=d2psidtau2(i-39,del2,tau2)
        derdr_helm=derdr_helm+(n(i)*del2*((d2bidtt*psi1)+ &
        (2.d0*dbidt*dsidt)+((capdel1**bco2(i-39))*d2sidtt)))
      enddo
      
      dpdtt = derdr_helm
      
end subroutine dphirdtautau

! ************************************************************************** !

subroutine dphirddeldtau(dpddt,del2,tau2)

      implicit none
      integer :: i
      real*8 :: del2,tau2,derdr_helm
      real*8 :: dpddt
      real*8 :: psi1,capdel1,dsidt,dsiddt,dbidd,dbidt,dsidd,d2biddt

!     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)*d(i)* &
        (tau2**(ti(i)-1))*(del2**(d(i)-1)))
      enddo
      
      do i = 8, 34
        derdr_helm=derdr_helm+(n(i)*ti(i)*(del2**(d(i)-1))* &
        (tau2**(ti(i)-1))*exp(-del2**c(i))* &
        (d(i)-(c(i)*(del2**c(i)))))
      enddo
      
      do i = 35, 39
        derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))* &
        exp((-alpha(i-34)*((del2-epsilon(i-34))**2.0d0))-(beta(i-34)* &
        ((tau2-gamma(i-34))**2.0d0)))* &
        ((ti(i)/tau2)-(2.0d0*beta(i-34)*(tau2-gamma(i-34))))* &
        ((d(i)/del2)-(2.d0*alpha(i-34)*(del2-epsilon(i-34)))))
      enddo
      
      do i = 40, 42
        psi1=psi(i-39,del2,tau2)
        capdel1=capdel(i-39,del2,tau2)
        dsidt=dpsidtau(i-39,del2,tau2)
        dsiddt=d2psiddeltau(i-39,del2,tau2)
        dbidd=ddelbiddel(i-39,del2,tau2)
        dbidt=ddelbidtau(i-39,del2,tau2)
        dsidd=dpsiddel(i-39,del2,tau2)
        d2biddt=d2delbiddeltau(i-39,del2,tau2)
        derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i-39))*(dsidt+ &
        (del2*dsiddt)))+(del2*dbidd*dsidt)+(dbidt* &
        (psi1+(del2*dsidd)))+(d2biddt*del2*psi1)))
      enddo
      
      dpddt = derdr_helm
      
end subroutine dphirddeldtau

! ************************************************************************** !

function dpsiddel(i,del2,tau2)
      implicit none
      real*8 :: dpsiddel
      real*8 del2,tau2,psi1
      integer i

      psi1=psi(i,del2,tau2)
      dpsiddel = -2.d0*capc(i)*(del2-1.d0)*psi1

      end function dpsiddel

! ************************************************************************** !

      function d2psiddel2(i,del2,tau2)
      implicit none
      real*8 :: d2psiddel2
      real*8 :: del2,tau2,psi1
      integer :: i

      psi1=psi(i,del2,tau2)
      d2psiddel2 = (2.d0*capc(i)*(del2-1.d0)**2-1.d0)*2.d0*capc(i)*psi1

      end function d2psiddel2

! ************************************************************************** !

      function dpsidtau(i,del2,tau2)
      implicit none
      real*8 :: dpsidtau
      real*8 :: del2,tau2,psi1
      integer :: i

      psi1=psi(i,del2,tau2)
      dpsidtau = -2.d0*capd(i)*(tau2-1.d0)*psi1

      end function dpsidtau

! ************************************************************************** !

      function d2psidtau2(i,del2,tau2)
      implicit none
      real*8 :: d2psidtau2
      real*8 :: del2,tau2,psi1
      integer :: i

      psi1=psi(i,del2,tau2)
      d2psidtau2 = ((2.d0*capd(i)*(tau2-1.d0)**2)-1.d0)*2.d0* &
      capd(i)*psi1

      end function d2psidtau2

! ************************************************************************** !

      function d2psiddeltau(i,del2,tau2)
      implicit none
      real*8 :: d2psiddeltau
      real*8 :: del2,tau2,psi1
      integer :: i

      psi1=psi(i,del2,tau2)
      d2psiddeltau = 4.d0*capc(i)*capd(i)*(del2-1.d0)*(tau2-1.d0)*psi1

      end function d2psiddeltau

! ************************************************************************** !

      function theta(i,del2,tau2)
      implicit none
      real*8 :: theta
      real*8 :: del2,tau2
      integer :: i

      theta=1.d0-tau2+capa(i)*((del2-1.d0)**2)**(1.d0/(2.d0*beta(i+5)))

      end function theta

! ************************************************************************** !

      function capdel(i,del2,tau2)
      implicit none
      real*8 :: capdel
      real*8 :: theta1,del2,tau2
      integer :: i

      theta1=theta(i,del2,tau2)
      capdel=theta1*theta1+capb(i)*(((del2-1.d0)**2)**aco2(i))
      end function capdel

! ************************************************************************** !

      function dcapdelddel(i,del2,tau2)
      implicit none
      real*8 :: dcapdelddel
      real*8 :: del2,tau2,theta1
      integer :: i

      theta1=theta(i,del2,tau2)
      dcapdelddel=(del2-1.d0)*((capa(i)*theta1*(2.d0/beta(i+5)) &
      *((del2-1.d0)**2 &
      )**((1/(2.d0*beta(i+5)))-1.d0))+(2.d0*capb(i)*aco2(i)* &
      (((del2-1.d0)**2)**(aco2(i)-1.d0))))

      end function dcapdelddel

! ************************************************************************** !

      function d2capdelddel2(i,del2,tau2)
      implicit none
      real*8 :: d2capdelddel2
      real*8 :: tmp1,del2,tau2,theta1
      real*8 :: ddd
      integer :: i

      theta1=theta(i,del2,tau2)
      tmp1=4.d0*capb(i)*aco2(i)*(aco2(i)-1.d0)*(((del2-1.d0) &
      **2)**(aco2(i)-2.d0))
 
      tmp1=tmp1+(2.d0*capa(i)*capa(i)*((1.d0/beta(i))**2)* &
      ((((del2-1.d0)**2)**((1.d0/(2.d0*beta(i)))-1.d0))**2))

      tmp1=tmp1+capa(i)*theta1*(4.d0/beta(i))*((1.d0/(2.d0*beta(i))) &
      -1.d0)*(((del2-1.d0)**2)**((1.d0/(2.d0*beta(i)))-2.d0))

      tmp1=tmp1*((del2-1.d0)**2.d0)

      ddd=dcapdelddel(i,del2,tau2)
      d2capdelddel2=tmp1+((1.d0/(del2-1.d0))*ddd)

      end function d2capdelddel2

! ************************************************************************** !

      function ddelbiddel(i,del2,tau2)
      implicit none
      real*8 :: ddelbiddel
      integer :: i
      real*8 :: del2,tau2,capdel1,ddd
      
      capdel1=capdel(i,del2,tau2)
      ddd=dcapdelddel(i,del2,tau2)
      ddelbiddel=bco2(i)*(capdel1**(bco2(i)-1.d0))*ddd

      end function ddelbiddel

! ************************************************************************** !

      function d2delbiddel2(i,del2,tau2)
      implicit none
      real*8 :: d2delbiddel2
      integer :: i
      real*8 :: del2,tau2,ddd1,ddd2,capdel1

      ddd1=dcapdelddel(i,del2,tau2)
      ddd2=d2capdelddel2(i,del2,tau2)
      capdel1=capdel(i,del2,tau2)
      d2delbiddel2=bco2(i)*ddd2*capdel1**(bco2(i)-1.d0) &
      +(bco2(i)-1.d0)*(capdel1**(bco2(i)-2.d0))*ddd1**2

      end function d2delbiddel2

! ************************************************************************** !

      function ddelbidtau(i,del2,tau2)
      implicit none
      real*8 :: ddelbidtau
      integer :: i
      real*8 :: del2,tau2,theta1,capdel1

      theta1=theta(i,del2,tau2)
      capdel1=capdel(i,del2,tau2)

      ddelbidtau=-2.d0*theta1*bco2(i)*capdel1**(bco2(i)-1.d0)

      end function ddelbidtau

! ************************************************************************** !

      function d2delbidtau2(i,del2,tau2)
      implicit none
      real*8 :: d2delbidtau2
      real*8 :: del2,tau2
      integer :: i
      real*8 :: capdel1,theta1

      capdel1=capdel(i,del2,tau2)
      theta1=theta(i,del2,tau2)
      d2delbidtau2=(2.d0*bco2(i)*(capdel1**(bco2(i)-1.d0)))+(4.d0* &
      (theta1**2)*bco2(i)*(bco2(i)-1.d0)*(capdel1**(bco2(i)-2.d0)))

      end function d2delbidtau2

! ************************************************************************** !

      function d2delbiddeltau(i,del2,tau2)
      implicit none
      real*8 :: d2delbiddeltau
      real*8 :: tmp3,del2,tau2
      integer :: i
      real*8 :: capdel1,ddd,theta1

      capdel1=capdel(i,del2,tau2)
      theta1=theta(i,del2,tau2)
      ddd=dcapdelddel(i,del2,tau2)
      tmp3=-capa(i)*bco2(i)*(2.d0/beta(i))*(capdel1**(bco2(i)-1.d0))* &
      (del2-1.d0)*(((del2-1.d0)**2)**((1.d0/(2.d0*beta(i)))-1.d0))

      tmp3=tmp3-(2.d0*theta1*bco2(i)*(bco2(i)-1.d0)*(capdel1** &
      (bco2(i)-2.d0))*ddd)
    
      d2delbiddeltau=tmp3

end function d2delbiddeltau

! ************************************************************************** !

subroutine vappr(tm,ps,dertp,derpt,ifl1)

      implicit none
      
!     co2 vapor pressure curve (Span & Wagner, 1996, p. 1524, eq. (3.13)

      integer :: j,maxit
      real*8 :: tm,ps,nu,dertp,derpt,ps1,uguess,lguess,xacc
      real*8 :: a1,a2,a3,a4
      integer :: ifl1
      real*8 :: dnu,ps2,f,df

      parameter(maxit=1000)
      parameter(xacc=1.d-7)
      a1 = -7.0602087d0
      a2 =  1.9391218d0
      a3 = -1.6463597d0
      a4 = -3.2995634
      if (ifl1 .eq. 11) then
        nu = 1.d0-(tm/tc)

        ps1 = a1*nu+a2*nu**1.5d0+a3*nu**2+a4*nu**4.d0
        ps2 = a1+1.5d0*a2*nu**0.5d0+2.0d0*a3*nu+4.d0*a4*nu**3

        ps = ps1*tc/tm

        ps = exp(ps)*pc

        dertp = (ps/tm)*(-ps2-((ps1*tc)/tm))
      else
        lguess = 50.
        uguess = 400.
        tm = (uguess-lguess)/2.d0
        nu = 1.d0-(tm/tc)            
        do j = 1, maxit
          ps1 = a1*nu+(a2*(nu**1.5d0))+(a3*(nu**2))+ &
                (a4*(nu**4))
          ps2 = a1+(1.5d0*a2*(nu**0.5d0))+(2.0d0*a3*nu)+ &
                (4.d0*a4*(nu**3))
          f = log(ps/pc)-(ps1/(1.d0-nu))
          df=-(ps1+(ps2*(1.d0-nu)))/((1.d0-nu)**2)
          dnu=f/df
          nu=nu-dnu
          if(abs(dnu).lt.xacc) goto 10
        enddo
 10     continue
        tm = tc*(1.d0-nu)
        ps1 = a1*nu+(a2*(nu**1.5d0))+(a3*(nu**2))+ &
             (a4*(nu**4.d0))
        ps2 = a1+(1.5d0*a2*(nu**0.5d0))+(2.0d0*a3*nu)+ &
             (4.d0*a4*(nu**3))
        derpt = (ps/tm)*(-ps2-((ps1*tc)/tm))
        derpt=1/derpt
      endif
      
end subroutine vappr

! ************************************************************************** !

subroutine viscosity(p,t,rho,drhodp,drhodt,mu,dmudt,dmudp)

! Fenghour, A., W. A. Wakeham, and V. Vesovic, 
! The viscosity of carbon dioxide, 
! J. Phys. Chem. Ref. Data, 27(1), 31−44, 1998.

      implicit none
      real*8 :: p, t, rho
      real*8 :: xsection1, xsection2, drhodt, drhodp,dmudp,dmudt
      real*8 :: dtxsection,dt_zerodenmu,dp_zerodenmu,drho_excessmu
      real*8 :: dp_excessmu, dt_excessmu
      real*8 :: t_star, xsection, zeroden_mu, excess_mu, mu
      real*8 :: lnstr,r2,r4,r5,r6,r7,r8

!     zero density viscosity
!     av(0) = 0.235156d0
!     av(1) = -.491266d0
!     av(2) = 5.211155d-2
!     av(3) = 5.347906d-2
!     av(4) = -1.537102d-2

      t_star = t/251.196
      
      lnstr = log(t_star)

      xsection = av(0) + (av(1) + (av(2) + (av(3) + &
          av(4)*lnstr)*lnstr)*lnstr)*lnstr
      
      xsection1 = xsection
      xsection = exp(xsection)

      zeroden_mu = 1.00697d0*(t**0.5d0)/xsection
      
      r2 = rho*rho
      r4 = r2*r2
      r5 = rho*r4
      r6 = rho*r5
      r7 = rho*r6
      r8 = r2*r6

      excess_mu = 0.4071119d-2*rho+0.7198037d-4*r2+ &
          0.2411697d-16*r6/t_star**3+ &
          0.2971072d-22*r8+ &
          (-0.1627888d-22*r8/t_star)

      mu = zeroden_mu+excess_mu

      xsection2 = av(1) + (2.d0*av(2) + &
          (3.d0*av(3) + 4.d0*av(4)*lnstr)*lnstr)*lnstr

      dtxsection = xsection2/(251.196d0*t_star)
      dt_zerodenmu =  exp(xsection1)*dtxsection

      dp_zerodenmu = 0.d0

      drho_excessmu = 0.4071119d-2+2.d0*0.7198037d-4*rho+ &
          6.d0*0.2411697d-16*r5/(t_star**3)+ &
          8.d0*0.2971072d-22*r7+ &
          8.d0*(-0.1627888d-22)*r7/t_star

      dp_excessmu = drho_excessmu

      dt_excessmu = 0.4071119d-2*drhodt+ &
          2.d0*0.7198037d-4*rho*drhodt+ &
          (0.2411697d-16*r5*((6.d0*t_star*drhodt) &
          -(3.d0*rho/251.196d0))/ &
          (t_star**4))+ &
          (8.d0*0.2971072d-22*r7*drhodt)+ &
          (-0.1627888d-22*r7*((8.d0*t_star*drhodt)- &
          (rho/251.196d0))/(t_star**2))

      dmudt = dt_excessmu + dt_zerodenmu

      dmudp = dp_zerodenmu + (dp_excessmu*drhodp)

end subroutine viscosity

! ************************************************************************** !

subroutine dissco2(p,t,mco2,fg,mol)

      implicit none
      real*8 p, t, mco2, fg, mol
      real*8 pc_h2o, tc_h2o, t1,ph2o,liq_cp, lambdaco2_na
      real*8 tauco2_na_cl, rhs

!     calculate solubility of CO2 based on Duan model

!     partial pressure of water based on empirical model provided
!     in Duan paper, Eqn. B1

!     convert pressure in bar and temperature in K

      p = p*10.d0
      pc_h2o = 220.85d0
      tc_h2o = 647.29d0

      t1 = (t/tc_h2o)-1.d0

      ph2o = (pc_h2o*t/tc_h2o)*(1.d0+(-38.640844d0*((-t1)**1.9d0)) &
          +(5.894842d0*t1)+(59.876516*t1*t1)+(26.654627*(t1**3.d0)) &
          +(10.637097*(t1**4.d0)))

      liq_cp = 28.9447706d0 + (-0.0354581768d0*t) + (-4770.67077d0/t) &
      +(1.02782768e-5*t*t) + (33.8126098/(630-t)) + (9.0403714d-3*p) &
      +(-1.14934031d-3*p*log(t)) + (-0.307405726*p/t) &
      + (-0.0907301486*p/(630.d0-t)) &
      + (9.32713393d-4*p*p/((630.d0-t)**2.d0))

      lambdaco2_na = -0.411370585d0 + (6.07632013d-4*t) &
      + (97.5347708d0/t) + (-0.0237622469d0*p/t) &
      + (0.0170656236*p/(630.d0-t)) + (1.41335834d-5*t*(log(p)))

      tauco2_na_cl = 3.36389723d-4 + (-1.9829898d-5*t) &
      + (2.1222083d-3*p/t) + (-5.24873303d-3*p/(630.d0-t))

      rhs = liq_cp - log(fg) + (2.d0*lambdaco2_na*mol) &
      + (tauco2_na_cl*mol*mol)

      mco2 = (p-ph2o)/exp(rhs)

      p = p*0.1d0

end subroutine dissco2

end module span_wagner_module
