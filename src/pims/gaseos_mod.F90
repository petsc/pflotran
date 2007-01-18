module Gas_Eos_Module
implicit none


  public
 
 
contains

  subroutine ideal_gaseos_noderiv(p,tc,energyscale,d,h,u)
    
    real*8,intent(in):: p,tc,energyscale
   real*8, intent(out):: d,h,u

   
    real*8, parameter:: Rg=8.31415D0 
    real*8, parameter:: Cpg=Rg*2.5 !maybe
    real*8  t
    
    t = tc + 273.15
    d= p / t / Rg /1D3
    h= Cpg * t * energyscale*1D3
    u= (Cpg- Rg) * t * energyscale*1D3

   
  end subroutine ideal_gaseos_noderiv

 
  subroutine ideal_gaseos(p,tc,energyscale,d,d_p,d_t,h,h_p,h_t,u,u_p,u_t)
    
    real*8, intent(in):: p,tc,energyscale
    real*8, intent(out):: d,d_p,d_t,h,h_p,h_t,u,u_p,u_t

    real*8, parameter:: Rg=8.31415 
    real*8, parameter:: Cpg=Rg*2.5 ! maybe
    real*8  t

    t = tc + 273.15
    d= p / t / Rg/1D3
    h= Cpg * t * energyscale*1D3
    u= (Cpg- Rg) * t * energyscale*1D3

    d_p=  d / p
    d_t=- d / t
    h_p=0.
    h_t=Cpg * energyscale*1D3
    u_p=0.
    u_t= (Cpg- Rg) * energyscale*1D3

!print *,'ideal gas ',energyscale,h,h_p,h_t,u,u_p,u_t
  end subroutine ideal_gaseos






!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!c   REFERENCES
!
!c     THIS ROUTINE IS LARGELY ADAPTED FROM THE TOUGH CODE.
!
!c     this routine computes the viscosity of vapor-air mixtures.
!c     it uses a modified version of a formulation based on kinetic
!c     gas theory, as given by j.o. hirschfelder, c.f. curtiss, and
!c     r.b. bird, molecular theory of gases and liquids, john wiley
!c     & sons, 1954, pp. 528-530.
!c
!c     the modification made to the hirschfelder et al. expressions is
!c     that for vapor viscosity accurate (empirical) values are used,
!c     rather than the first order expression of kinetic theory.
!c
!c     the formulation matches experimental data on viscosities of
!c     vapor-air mixtures in the temperature range from 100 to 150
!c     deg. c, for all compositions, to better than 4%.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine visgas_noderiv(t,pa,p,ds,visg)
      real*8  t,pa,p,ds,visg
      real*8  fmwh2o,fmwa,fair,fwat,cair,cwat

      data  fmwh2o,    fmwa,  fair,   fwat,    cair,  cwat &
           /18.0153d0, 28.96d0, 97.d0, 363.d0, 3.617d0, 2.655d0/
 
      real*8 fmix,cmix,d,xga,xg1,tk,trd1,trd3,ome1,ome3,ard,fmw3,vis1, &
             v1,vs,vis2,vis3,z1,g,h,e,z2,z3


!c======================================================================

      fmix = sqrt (fair*fwat)
      cmix = (cair+cwat)*0.5d0

!      do k = 1,nb
 !       if (iphas(k).eq.2 .or. iphas(k).eq.0) then

          d   = ds *fmwa       
          xga =0.5 ! pa /p ! for debug, set x constant
          xg1 = 1.D0 - xga
          tk  = t +273.15d0

          trd1 = tk/fair
          trd3 = tk/fmix
          ome1 = (1.188d0-0.051d0*trd1)/trd1
          ome3 = (1.480d0-0.412d0*log(trd3))/trd3
          ard  = 1.095d0/trd3
          fmw3 = 2.d0*fmwa*fmwh2o/(fmwa+fmwh2o)
          vis1 = 266.93d-7*sqrt(fmwa*trd1*fair)/(cair*cair*ome1*trd1)
 
          v1 = .407d0*t +80.4d0
          if (t .le.350.d0) then
            vs = 1.d-7*(v1-d*(1858.d0-5.9d0*t )*1.d-3)
          else
!             if (t .gt.350.d0) 
!cpcl .      vs = 1.d-7*(v1 + 0.353d0*d + 676.5d-6*d**2 + 102.1d-9*d**3)
           vs = 1.d-7*(v1 + (0.353d0 + (676.5d-6 + 102.1d-9*d)*d)*d)
          endif

          vis2 = 10.d0*vs
          vis3 = 266.93d-7*sqrt(fmw3*trd3*fmix)/(cmix*cmix*ome3*trd3)
          z1   = xga*xga/vis1+2.d0*xg1*xga/vis3+xg1*xg1/vis2
          g    = xga*xga*fmwa/fmwh2o
          h    = xg1*xg1*fmwh2o/fmwa
          e    = (2.d0*xga*xg1*fmwa*fmwh2o/fmw3**2)*vis3/(vis1*vis2)
          z2   = 0.6d0*ard*(g/vis1+e+h/vis2)
          z3   = 0.6d0*ard*(g+e*(vis1+vis2)-2.d0*xga*xg1+h)
          visg  = (1.d0+z3)/(z1+z2)*.1d0
           
      return
    end subroutine visgas_noderiv



  subroutine Henry_air_noderiv(p,tc,ps,Henry)
! Calculate Henry Coefficient for N2
! t in K
! Henry have the same unit as p and ps, then make it dimensionless by
! devide it with p

    implicit none
    real*8,intent(in) ::  p,tc,ps
    real*8,intent(out)::  Henry

    real*8  Tr,tao,tmp,t
    real*8, parameter :: a=-9.67578, b=4.72162, c=11.70585
    real*8, parameter :: Tcl=647.096 ! H2O critical temp(K) from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + B * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

   return 
  end subroutine Henry_air_noderiv


 subroutine Henry_air(p,tc,ps,ps_p,ps_t,Henry,Henry_p,Henry_t)
   implicit none
    real*8,intent(in) ::  p,tc,ps,ps_p,ps_t
    real*8,intent(out)::  Henry,Henry_p,Henry_t
! note t/K, p/Pa, Henry/Pa 

    real*8  Tr,tao,tmp,t
    real*8, parameter :: a=-9.67578, b=4.72162, c=11.70585
    real*8, parameter :: Tcl=647.096 ! H2O critical temp from IAPWS(1995b)

    t=tc+273.15D0
    Tr=t/Tcl
    tao=1.D0-Tr
    tmp= a/Tr + b * tao**0.355/Tr + c * (Tr**(-0.41)) * exp(tao)
    Henry=exp(tmp)*ps

    tmp =((-a/Tr+b*(-0.355*tao**(-0.645)-tao**0.355/Tr))/Tr - &
         c*exp(tao)*(tao**(-.41))*(0.41/Tr-1.))/Tcl
    Henry_t=Henry*(tmp +ps_t/ps)
    Henry_p=ps_p*Henry/ps

  
   return 
 end subroutine Henry_air

end module Gas_Eos_Module
