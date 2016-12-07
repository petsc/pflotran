! Modulex will act as a shield to seprate equation solving part
!  from the parameter evaluation part


! phase index:: 1: liquid,water rich, 2: Gas phase, 3: super critical CO2
! component index::1 h2o, 2 air, 3 co2   
! diffusion coefficient index:: (iphase-1)*num_comp+icomp
   

module mixture_module
#include "finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

 implicit none

 public
  
 PetscReal, private, parameter:: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                              fmwco2 = 44.0098D0
 PetscReal, private, parameter::eps=1.D-6

 contains

! ************************************************************************** !

subroutine mixture_eos_noderiv (p,t,xga,sg,energyscale,num_phase,num_spec,&
                    num_pricomp,ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,sat_pressure,den,avgmw,h, &
                    u,diff,hen,xphi,pc,kvr,ierr,itable)
  ! 
  ! subroutines to calculate the properties of mixture
  ! will:: call other EOS mod to obtain pure fluid properties
  ! apply mixing rules
  ! 
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
    use Water_EOS_module
    use Gas_EOS_module  
    use pckr_module
    use co2eos_module
    use co2_span_wagner_module


    implicit none

    PetscReal :: p,t,xga,sg,energyscale
    PetscInt :: ipckrtype !, ithrmtype
    PetscInt :: num_phase,num_spec,num_pricomp
    PetscReal :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
    PetscReal :: sat_pressure
    PetscReal :: den(num_phase),h(num_phase),u(num_phase), &
              avgmw(num_phase),pc(num_phase),kvr(num_phase)
    PetscReal :: diff(1:num_phase*num_spec),hen(1:num_phase*num_spec)
    PetscReal :: tmp

    PetscErrorCode :: ierr
    PetscInt :: iitable
    PetscInt, optional :: itable
    PetscReal :: tk, kr(num_phase),xla,xgw,xlw
    PetscReal :: dstea,dsteamol,dstea_p,dstea_t,hstea,hstea_p,hstea_t
    PetscReal :: dg,hg,visg,ug
    PetscReal :: pw,dw_kg, dw_mol,hw,visl
    PetscReal :: dif(1:num_spec),henry,pa
    PetscReal :: co2_poyn,xphi
    PetscReal :: dddt,dddp,fg,dfgdp,dfgdt,eng,dhdt,dhdp,dvdt,dvdp
    PetscInt :: iflag

    ierr=0; iitable=0
    if(present(itable)) iitable = itable
    h=0.D0; u=0.D0; avgmw=0.D0; pc=0.D0; kvr=0.D0
    diff=0.D0; hen=0.D0
 
  
    
    tk = t + 273.15D0

    call PSAT(t, sat_pressure, ierr)




!***************  Liquid phase properties **************************

   ! pure water
    pw=p
    if(num_phase>=2)then
      call pflow_pckr_noderiv_org(ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
                              pckr_m,pckr_pcmax,sg,pc,kr,pckr_betac,pckr_pwr)
      pw=p !-pc(1)
      if((pw<0.D0).and.(sg>(1.D0-eps)))then
        pw=p
        pc(1)=0.D0
      endif
    end if

    call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
   !    call VISW(t,pw,sat_pressure,visl,tmp,tmp2,ierr)
   ! call VISW_FLO(t,dw_mol,visl)
    call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
  !  print *,'visw  ',visl,tmp
  ! dif= 1.D-7 !effective diffusion coeff in liquid phase: check pcl
    dif= 1.D-9 ! m^2/s @ 25C

    !apply mixing rules
    ! should be careful with the consistance of mixing rules
    den(1) = dw_mol  !*c+(1-c)* 
    h(1) = hw 
    u(1) = h(1) - pw/den(1)* energyscale
    diff(1:num_spec) = dif
    kvr(1) = kr(1)/visl
    xlw = 1.D0
    avgmw(1) = xlw*fmwh2o+(1.d0-xlw)*fmwco2

!*****************Gas phase properties 2*************************
  if((num_phase>=2)) then
        
 
!********************************************************************* 
    
    pa=p*xga
    if(pa>=5d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      dg=den(2)*fmwco2
      iflag = 1
      call co2_span_wagner(p*xga*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag,iitable)
      dg= dg / fmwco2
      fg= fg * 1.D6
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
      fg=pa
    endif
!    call duanco2(t,p*xga/1D5,dg,fugco2,co2_phi)
!    dg= dg * 1.D3 /fmwco2
!    call ENTHALPY(T+273.15D0,1.D-3/dg ,1.D0/co2_phi, hg)
!    hg=hg*1.D-3

    xphi = fg/(p*xga)

    call Henry_CO2_noderiv(xla,tmp,t,p*xga,xphi,henry,co2_poyn)

!   call Henry_air_noderiv(p,t,sat_pressure,Henry)

!   if(sg<1.D0-eps)then
  
      hen(2)=p/henry  ! rkh [bars] Henry constant
!     hen(num_spec+1:num_spec*2)=1.0
!     hen(2)=1.D-4    !p/Henry
      hen(num_spec+1:num_spec*2)=0.0

!   else
  !   hen(2)=1.
  !   hen(num_spec+1:num_spec*2)=1.0
  ! endif
   

    pa = p * xga
    xgw = 1.d0-xga
    xla = hen(2)*xga
    xlw = 1.d0-xla



!*************** Checking EOS validation*******

!    if((p-pa)>(1.05D0*sat_pressure) .and.sg>(1.D0-eps) )then
     ! ierr=-1;
 !     print *,'Varible out of Reasonable Region G->L', p,p-pa,sat_pressure,(p-pa)/sat_pressure
     ! return
  !    elseif((p-pa)<0.95*sat_pressure .and.sg<eps)then
      !   ierr=-1;
   !    print *,'Varible out of Reasonable Region L->G', p,p-pa,sat_pressure,(p-pa)/sat_pressure
      ! return
  !  end if

!****************End Checking *****************
  
  !     call ideal_gaseos_noderiv(p,t,energyscale,dg,hg,ug)

!****   Steam Properties ***********************
!  if((p*xgw)<=sat_pressure)then 

    if(xgw>eps)then
       call steameos(t,p,pa,dstea,dsteamol,dstea_p,dstea_t,&
            hstea,hstea_p,hstea_t,energyscale,ierr) 
!print *,'steameos', t,p,pa,p-pa,sat_pressure, dsteamol,hstea,dg,hg,dsteamol/xgw
       dsteamol=dsteamol
       hstea=hstea!/xgw
!       tmp=dsteamol
    else
        dsteamol=dg /xga*xgw
        hstea=hg 
    end if


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


 den(2)= dg + dsteamol
   
!   call visgas_noderiv(t,pa,p,den(2),visg)
!call visgas_noderiv(t,pa,p,den(2),visg)
! call visco2(t,dg ,visg)
!    visg=8.d-6
    kvr(2)=kr(2)/visg

  
!    den(2)=1.D0/( xga/dg + xgw/dsteamol)  !*c 
  !  den(2)=1.D0/( xga/dg + 1.D0/dsteamol)
    avgmw(1)=xlw * fmwh2o + xla * fmwco2
    avgmw(2)=xga*fmwco2+xgw*fmwh2o
    h(2)=  hg *xga + hstea*xgw 
 !   h(2)= ( hg*xga  + hstea*xgw ) 
 !   h(2)= ( hg *xga + hstea ) 
    u(2)=  h(2)-p/den(2) * energyscale
    pc(2)=0
   !  print *,'gas phase nonder h::',t,p,h(2),u(2),hg,hstea
    diff(num_spec+1:num_spec*2)=1.D-5 ! dirty part
    
 !if(t>=100.D0) print *, p,t,xga,xgw,dg,dsteamol,hg, hstea, h(2)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! remember: henry coeff is arranged in the vector as 
!    phase          1(water)          2(gas)            3(co2)
!   species      1    2    3       1    2    3        1   2    3
!   values                        1.0  1.0  1.0
!_________________________________________________________________
 endif

  

 if(num_phase>=3)then 
    ! call
    ! call
    !  den(2)=dw_mol  !*c+(1-c)* 
    !    h(2)=  hw
    !    u(2)=  hw-p
    !    pc(2)=0
    !    kvr(2)=kr/vis
 endif
 return
end subroutine mixture_eos_noderiv

! ************************************************************************** !

subroutine mixture_eos(p,t,xga,sg,energyscale,num_phase,num_spec, &
    num_pricomp,ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,pckr_m, &
    pckr_pcmax,pckr_betac,pckr_pwr,sat_pressure,den,den_p,den_t,den_c,den_s, &
    avgmw,avgmw_c, &
    h, h_p, h_t, h_c, h_s, u, u_p, u_t, u_c, u_s, diff, diff_p, &
    diff_t, diff_c, diff_s, hen, hen_p, hen_t, hen_c, &
    hen_s, xphi, pc, pc_p, pc_t, pc_c, pc_s, kvr, kvr_p, kvr_t, kvr_c, kvr_s, &
    ierr,itable)


    use Water_EOS_module
    use Gas_EOS_module  
    use pckr_module
    use co2eos_module
    use co2_span_wagner_module 

    implicit none

! Notice: kvrl=krl/visl  kvrg=krg/visg
!  since the conservation equations do not refer to kr or vis  individually
    PetscReal  :: p,t,xga,sg,energyscale
    PetscInt :: num_phase,num_spec,num_pricomp, ipckrtype
    PetscInt, optional :: itable
    PetscReal  :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
    PetscReal  :: sat_pressure
    PetscReal  :: den(num_phase),den_p(num_phase),den_t(num_phase), &
                              den_s(num_phase),den_c(num_phase*num_pricomp)
    PetscReal  :: avgmw(num_phase),avgmw_c(num_phase*num_pricomp)
    PetscReal  :: h(num_phase),h_p(num_phase),h_t(num_phase), &
                            h_s(num_phase),h_c(num_phase*num_pricomp)
    PetscReal  :: u(num_phase),u_p(num_phase),u_t(num_phase), &
                            u_s(num_phase),u_c(num_phase*num_pricomp)
    PetscReal  :: pc(1:num_phase),pc_p(1:num_phase),pc_t(1:num_phase), &
                             pc_s(1:num_phase),pc_c(1:num_phase*num_pricomp)
    PetscReal  :: kvr(1:num_phase),kvr_p(1:num_phase),kvr_t(1:num_phase), &
                              kvr_s(1:num_phase),kvr_c(1:num_phase*num_pricomp)
    PetscReal  :: diff(1:num_phase*num_spec),diff_p(1:num_phase*num_spec),&
               diff_t(1:num_phase*num_spec),diff_s(1:num_phase*num_spec),&
               diff_c(1:num_phase*num_spec*num_pricomp)
    PetscReal  ::  hen(1:num_phase*num_spec),hen_p(1:num_phase*num_spec),&
               hen_t(1:num_phase*num_spec),hen_s(1:num_phase*num_spec),&
               hen_c(1:num_phase*num_spec*num_pricomp)

    PetscInt :: iitable
    PetscErrorCode :: ierr
    PetscReal :: sat_pressure_t,sat_pressure_p,pw,dw_kg,dw_mol,dw_p,dw_t,hw,visl
    PetscReal :: hw_p,hw_t,dg,hg,ug
!   PetscReal :: dg_p,dg_t,hg_p,hg_t,ug_p,ug_t
    PetscReal :: dstea,dsteamol,dstea_p,dstea_t, hstea,hstea_p,hstea_t
    PetscReal :: kr(num_phase),kr_s(num_phase),visg,visw_p,visw_t,visg_t,&
               visg_p,visg_c 
    PetscReal :: ps,pa,henry,henry_p,henry_t,henry_c,henry_s
    PetscReal :: xla,xgw,xlw,xlco2,xmlco2
    PetscReal :: tmp,tmp2,tmp3,tmp0,tmp4,tmp5
    PetscReal :: co2_poyn,xphi
    PetscReal :: dddt,dddp,fg,dfgdp,dfgdt,eng,dhdt,dhdp,dvdt,dvdp,tmpdg
    PetscInt :: iflag

    ierr=0; iitable=0
    if(present(itable)) iitable =itable
    den_p=0.D0; den_t=0.D0; den_c=0.D0; den_s=0.D0
    avgmw=0.D0; avgmw_c=0.D0 
    h=0.D0; h_p=0.D0; h_t=0.D0; h_c=0.D0; h_s=0.D0 
    u=0.D0; u_p=0.D0; u_t=0.D0; u_c=0.D0; u_s=0.D0
    diff=0.D0; diff_p=0.D0; diff_t=0.D0; diff_c=0.D0; diff_s=0.D0
    henry=0.D0; henry_p=0.D0; henry_t=0.D0; henry_c=0.D0; henry_s=0.D0
    pc=0.D0;  pc_p=0.D0;  pc_t=0.D0; pc_c=0.D0; pc_s=0.D0 
    kvr=0.D0; kvr_p=0.D0; kvr_t=0.D0; kvr_c=0.D0; kvr_s=0.D0
   
  
   

    !call PSAT(t, sat_pressure, ierr)
    !ps=sat_pressure
    call PSAT1(t,sat_pressure,sat_pressure_t,ierr)
    ps=sat_pressure
    sat_pressure_t =1.D0/sat_pressure_t !in PSAT1 it is dT/dPsat
    sat_pressure_p=0.D0

!***************  Liquid phase properties **************************
    pw=p
    if(num_phase>=2)then
       call pflow_pckr(ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
          pckr_m,pckr_pcmax,sg,pc,pc_s,kr,kr_s,pckr_betac,pckr_pwr) 
    !   print *,'kr', kr,kr_s
       pw=p!-pc(1)
       if((pw<0.D0).and.(sg>(1.D0-eps)))then
          pw=p
          pc(1)=0.D0
          pc_s(1)=0.D0
       endif
    endif
   ! pure water
    call wateos(t,pw,dw_kg,dw_mol,dw_p,dw_t,hw,hw_p,hw_t,energyscale,ierr)
    call VISW(t,pw,sat_pressure,visl,visw_t,visw_p,ierr)
  ! print *,'vis w:: ',pw,visl
    
 !  once called, gives pc, kr for both phases 
   
   
   
  
    !apply mixing rules
    ! should be careful with the consistance of mixing r

    den(1)=dw_mol ! should be changed to the same as gas phase
    den_p(1)= dw_p
    den_t(1)= dw_t
    den_c(1)= 0.D0
    den_s(1)= 0.D0  !-dw_p*pc_s(1)

    xlw=1.D0
    avgmw(1)=xlw*fmwh2o+(1.D0-xlw)*fmwco2
    avgmw_c(1)=0.D0 

    h(1)  = hw
    h_p(1)=hw_p
    h_t(1)=hw_t
    h_s(1)=0.D0 !-hw_p*pc_s(1)
  
    tmp=den(1)*den(1)
    u(1)=  h(1) - pw/den(1) * energyscale 
    u_p(1)=hw_p-(den(1)-pw*den_p(1))/tmp * energyscale 
    u_t(1)=hw_t + pw/tmp*den_t(1) * energyscale 
    u_s(1)=(-pc_s(1)/den(1)+p/tmp*den_s(1)) * energyscale 
    u_c(1)=0.d0

    diff(1:num_spec)= 1.D-7
    diff_p(1:num_spec)=0.d0
    diff_t(1:num_spec)=0.d0
    diff_s(1:num_spec)=0.d0
    diff_c(1:num_spec*num_pricomp)=0.d0


    kvr(1)=kr(1)/visl
    kvr_p(1)=-kr(1)/visl/visl*visw_p
    kvr_t(1)=-kr(1)/visl/visl*visw_t
    kvr_s(1)=kr_s(1)/visl
    kvr_c(1)=0.d0

  ! print *,'kvr', kvr(1),kvr_p(1), kvr_t(1),kr(1),visl,visw_p,visw_t

!*****************Gas phase properties 2*************************
 if((num_phase>=2))then
    pa=p * (1.D0-xgw)
    !call  ideal_gaseos(p,t,energyscale,dg,dg_p,dg_t,hg,hg_p,hg_t,ug,ug_p,ug_t)

!********************************************************************* 


    !call Henry_air(p,t,sat_pressure,sat_pressure_p,sat_pressure_t, Henry,&
    !          Henry_p,Henry_t)
    
    !call duanco2(t,p*xga/1D5,dg,fugco2,co2_phi)
    !dg= dg * 1.D3 /fmwco2
    !call ENTHALPY(T+273.15D0,1.D-3/dg ,1.D0/co2_phi, hg)
    !hg=hg*1.D-3
 
    pa=p*xga
    if(pa>=5.d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      dg=den(2)*fmwco2
      iflag = 1
      call co2_span_wagner(p*xga*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag,iitable)
      tmpdg = dg 
      dg= dg / fmwco2
      fg= fg * 1.D6
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
      fg=pa
    endif

    xphi = fg/(p*xga)

    call Henry_CO2_noderiv(xla,tmp,t,p*xga,xphi,henry,co2_poyn)

    hen(2) = p/Henry
    hen_t(2) = -hen(2)/Henry*Henry_t
    hen_p(2) = (1.D0 - hen(2)*Henry_p)/Henry

!   print *,'mixture_eos: ',p,fg,p*xga,t,xla,xga,henry,den(2),xphi

    hen(1)=0.D0
    hen_t(1)=0.D0
    hen_p(1)=0.D0
    hen(num_spec+1:num_spec*2)=0.D0

!****   Steam Properties ***********************
    xgw=1.D0-xga; xla=hen(2)*xga; xlw=1.D0-xla; pa=p*xga


!*************** Checking EOS validation*******
!print *,'vapor p :  ', p,p-pa,sat_pressure
 !   if((p-pa)>(1.05D0*sat_pressure).and.(sg>eps) .and.(sg<(1.D0-eps)))then
      !ierr=-1;
  !   print *,'Varible out of Reasonable Region der1',sg,p, (p-pa),sat_pressure,(p-pa)/sat_pressure,eps
      !return
   !  elseif((p)<0.95*sat_pressure .and.sg<(1.D0-eps).and.sg>eps)then
!  print *,'------------Varible out of Reasonable Region der2',sg,p, (p-pa),sat_pressure,(p-pa)/sat_pressure
 !   end if

!****************End Checking *****************
  

   
    if(xgw>eps)then
      call steameos(t,p,pa,dstea,dsteamol,dstea_p,dstea_t,&
            hstea,hstea_p,hstea_t,energyscale,ierr) 
      dsteamol=dsteamol!/xgw
      hstea=hstea!/xgw
      tmp=dsteamol
    else
      dsteamol=dg /xga*xgw
      hstea=hg
    end if

    den(2) = dg + dsteamol
    h(2) = hg*xga + hstea *xgw 
    u(2) = h(2)-p/den(2) * energyscale
   ! call visgas_noderiv(t,pa,p,den(2),visg)

  ! conduct numerical evaluation of vapor properties ********************************
    ! dp  
    tmp=-1.D2
   !call  ideal_gaseos(p+tmp,t,energyscale,dg,dg_p,dg_t,hg,hg_p,hg_t,ug,ug_p,ug_t) 

    pa=(p+tmp)*xga
    
    if(pa>=5.d4)then
      dg=tmpdg
      iflag = 1
      call co2_span_wagner((p+tmp)*xga*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,tmp4,dvdt,dvdp,iflag,iitable)
      dg= dg / fmwco2
      fg= fg * 1.D6
      hg= hg * fmwco2
    else ! use ideal gas
      call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,tmp4)
      fg=pa
    endif
   
    call Henry_CO2_noderiv(xlco2,xmlco2,t,(p+tmp)*xga,fg/(p+tmp)/xga, &
    henry,co2_poyn)

    if(xgw>eps)then
      call steameos(t,p+tmp,(p+tmp)*xga,dstea,dsteamol,dstea_p,dstea_t,&
               hstea,hstea_p,hstea_t,energyscale,ierr)
    else
      dsteamol=dg/xga*xgw
      hstea=hg 
    end if

    tmp0 = dg + dsteamol !den
    tmp2 = hstea*xgw + hg *xga                    !h
    tmp3 = tmp2 - (p+tmp)/tmp0 * energyscale      !u
    tmp5 = (p+tmp)/henry
       
    den_p(2)= (tmp0-den(2))/tmp 
    hen_p(2)=(tmp5-hen(2))/tmp
    h_p(2) = (tmp2-h(2))/tmp
    u_p(2) = (tmp3 - u(2))/tmp
    visg_p = (tmp4-visg)/tmp

    !dt
    tmp=1.D-3
     ! call  ideal_gaseos(p,t+tmp,energyscale,dg,dg_p,dg_t,hg,hg_p,hg_t,ug,ug_p,ug_t)
    pa=p*xga
    if(pa>=5d4)then
      dg=tmpdg
      iflag = 1
      call co2_span_wagner(p*xga*1.D-6,t+tmp+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,tmp4,dvdt,dvdp,iflag,iitable)
      dg= dg / fmwco2
      fg= fg * 1.D6
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(pa,t+tmp,energyscale,dg,hg,ug)
      call visco2(t+tmp,dg*fmwco2,tmp4)
      fg=pa
    endif

    call Henry_CO2_noderiv(xlco2,xmlco2,t+tmp,p*xga,fg/p/xga,henry,co2_poyn)
    if(xgw>eps)then
      call steameos(t+tmp,p,pa,dstea,dsteamol,dstea_p,dstea_t,&
             hstea,hstea_p,hstea_t,energyscale,ierr)
    else
      dsteamol=dg/xga*xgw
      hstea=hg
    end if
    tmp0 = dg + dsteamol !den
    tmp2 = hstea *xgw+ hg * xga                !h
    tmp3 = tmp2 - p/tmp0 * energyscale         !u
    tmp5 = p/henry
       
    den_t(2)= (tmp0-den(2))/tmp 
    hen_t(2)=(tmp5-hen(2))/tmp
    h_t(2) = (tmp2 - h(2))/tmp
    u_t(2) = (tmp3 - u(2))/tmp
    visg_t = (tmp4-visg)/tmp

    !dx   
    tmp=-1.D-6
    if(xga<0.1D0) tmp = 1.D-5
    xgw=1.D0-(xga+tmp); xla=hen(2)*(xga+tmp) ; xlw=1.D0 - xla;pa=p * (xga+tmp)
       !call  ideal_gaseos(p,t,energyscale,dg,dg_p,dg_t,hg,hg_p,hg_t,ug,ug_p,ug_t)

  !    call duanco2(t,p*(xga+tmp)/1D5,dg,fugco2,co2_phi)
  ! dg= dg * 1.D3 /fmwco2
  ! call ENTHALPY(T+273.15D0,1.D-3/dg ,1.D0/co2_phi, hg)
  ! hg=hg*1.D-3
  ! call Henry_CO2_noderiv(xlco2,xmlco2,t,p*(xga+tmp),co2_phi,henry,co2_poyn)

    pa=p*(xga+tmp)
    if(pa>=5d4)then
      dg=tmpdg
      iflag = 1
      call co2_span_wagner(p*(tmp+xga)*1.D-6,t+273.15D0,dg,dddt,dddp, &
      fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,tmp4,dvdt,dvdp,iflag,iitable)
      dg= dg / fmwco2
      fg= fg * 1.D6
      hg= hg * fmwco2
    else      
      call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
      call visco2(t,dg*fmwco2,tmp4)
      fg=pa
    endif

    call Henry_CO2_noderiv(xlco2,xmlco2,t,p*(tmp+xga),fg/p/(xga+tmp), &
    henry,co2_poyn)

    if((1.D0-xga)>eps)then
      call steameos(t,p,pa,dstea,dsteamol,dstea_p,dstea_t,&
               hstea,hstea_p,hstea_t,energyscale,ierr)
    else
      dsteamol=dg/(xga+tmp)*xgw
      hstea=hg
    end if
    tmp0 =  dg + dsteamol   !den
    tmp2 = hstea *xgw+ hg *(xga+tmp)               !h
    tmp3 = tmp2 - p/tmp0 * energyscale             !u
    tmp5 = p/henry           !henry
       
    den_c(2)= (tmp0-den(2))/tmp
    hen_c(2)= (tmp5-hen(2))/tmp
    h_c(2) = (tmp2 - h(2))/tmp
    u_c(2) = (tmp3 - u(2))/tmp
    visg_c = (tmp4-visg)/tmp

    xgw=1.D0-xga; xla=hen(2)*xga ; xlw=1.D0 - xla;pa=p * xga
    
! print *,'gas phase der h::',t,p,xga,h(2),u(2),hg,hstea,dsteamol
!  dsteamol=dsteamol/xgw
!  dstea_t=dstea_t/xgw
!  hstea=hstea/xgw
!  hstea_t=hstea_t/xgw

!    dsteamol=dg
!    dstea_p=dg_p
!    dstea_t=dg_t
!    hstea=hg*10.0D0! + 35.0D0
!    hstea_p=hg_p*10.0D0
!    hstea_t=hg_t *10.0D0
!  Now just make the steam same as ideal gas   
!**********************************************************************


  
  ! visg=8.D-6
  ! visg_p=0.
  ! visg_t=0.
  ! visg_c=0.
  ! visg_s=0.
!print *,'Mixing:: ', xla,xga,xgw,dg,dsteamol,den(2)
 
!    Henry=Henry/p
!    den(2) =1.D0/( xga/dg + xgw/dsteamol)
  !  tmp=(den(2)/dg); tmp=tmp*tmp
  !  tmp2=(den(2)/dsteamol) ; tmp2=tmp2*tmp2
  !  den_p(2) = xga * tmp * dg_p + xgw * tmp2 *  dstea_p
  !  den_t(2)= xga * tmp * dg_t + xgw * tmp2 *  dstea_t
  !  den_c(2)=den(2)*den(2)*(1.D0/dsteamol-1.D0/dg)
!print *,'Mixing :: ', xla,xga,xgw,dg,den(2),den_p(2)

  !  h(2)=  hg *xga + hstea  * xgw
  !  h_p(2)=hg_p*xga + hstea_p *xgw  
  !  h_t(2)=hg_t*xga + hstea_t * xgw 
  !  h_c(2)= (hg *xga + hstea  * xgw- hg *(xga-1d-5) - hstea  * (xgw+1d-5))&
   !         /1d-5
!    if(h_c(2)<-1d3) h_c(2)=-1d3
  ! print *,'gas phase der h::',t,p,h(2),h_p(2),h_t(2), h_c(2),hg,hstea

  !  tmp=den(2)*den(2)
  !  u(2)=  h(2) - p/den(2) * energyscale
  !  u_p(2)= h_p(2) - (1.D0/den(2)-p/tmp*den_p(2)) * energyscale
  !  u_t(2)= h_t(2) + p/tmp*den_t(2) * energyscale
  !  u_c(2)= h_c(2) + p/tmp*den_c(2) * energyscale
   ! print *,'gas u der :: ',t,p, u(2),u_p(2),u_t(2),u_c(2),den_c(2),h_c(2)
 
    diff(1+num_spec : 2 * num_spec)= 1.D-5
    diff_p(1+num_spec : 2 * num_spec)=0.D0 
    diff_t(1+num_spec : 2 * num_spec)=0.D0
    diff_s(1+num_spec : 2 * num_spec)=0.D0
    diff_c(1+num_spec*num_pricomp : 2 * num_spec*num_pricomp)=0.D0
    
    kvr(2)=kr(2)/visg
    kvr_p(2)=-kvr(2)/visg * visg_p
    kvr_t(2)=-kvr(2)/visg * visg_t
    kvr_s(2)=kr_s(2)/visg
    kvr_c(2)=-kvr(2)/visg * visg_c
   

!   print *,1.D0-sg,kr(1),kr(2),pc(1),kr_s(1),kr_s(2),pc_s(2)

    avgmw(1) = xlw*fmwh2o + xla* fmwco2
    avgmw(2) = xgw* fmwh2o+xga*fmwco2
    avgmw_c(2) = fmwco2-fmwh2o
    avgmw_c(1) = Hen(2) * avgmw_c(2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! remember: henry coeff is arranged in the vector as 
!    phase          1(water)          2(gas)            3(co2)
!   species      1    2    3       1    2    3        1   2    3
!   values                        1.0  1.0  1.0
!_________________________________________________________________
 end if
! print *,'mix: ',p,t, xga, u(1), den(2),h(2),u(2),u_p(2),u_t(2),u_c(2)
end subroutine mixture_eos

end module mixture_module
