!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!c VERSION/REVISION HISTORY
! 
!c $Id: flpckr.f,v 1.1.1.1 2002/01/18 23:53:37 lichtner Exp $
!c $Log: flpckr.f,v $
!c Revision 1.1.1.1  2002/01/18 23:53:37  lichtner
!c Initial Entry
!c
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


module pckr_module
  public

  real*8, private, parameter:: pckr_sat_water_cut = 1.D0 - 5.D-7


  contains 

      subroutine pflow_pckr_noderiv_org(ipckrtype,pckr_swir,pckr_lambda, &
                 pckr_alpha,pckr_m,pckr_pcmax,sg,pc,kr,pckr_beta,pckr_pwr) 
      
      implicit none 


      integer ipckrtype
   !formation type, in pflow should be refered by grid%icap_loc
      real*8 :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax
      real*8 :: pckr_beta,pckr_pwr
      real*8 :: sg
      real*8 :: pc(1:2),kr(1:2)
       
      real*8 :: se,swir,sw0,lam,ala,um,un,upc,upc_s,kr_s, krg_s
      real*8 :: temp,ser,pcmax,sw
      real*8 :: uum,pckr_betac,betac,st
     
     
    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
      sw=1.D0-sg
     
      sw0=1.d0
      pcmax=pckr_pcmax
      swir=pckr_swir
      select case(ipckrtype)

      case(1) ! van Gennuchten
          ala=pckr_alpha
      !    swir=pckr_swir
          um=pckr_m
          un=1.D0/(1.D0-um)
          if(sw> pckr_sat_water_cut)then
            upc=0.D0; kr(1)=1.d0; kr(2)=0.d0;
          elseif(sw>(1.05D0*swir))then
            se=(sw-swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.d0
            !kr(2)=1.D0-kr(1)
      kr(2)=sqrt(1.D0 - se)*((1.D0-se**(1.D0/um))**um)**2.D0
       !    print *,'in pckr nond ',sw,se,upc,kr
          else  ! use linear extropolation
            se=(0.05D0*swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
            if(sw>swir)then
              kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.D0
        kr(2)=sqrt(1.D0 - se)*((1.D0-se**(1.D0/um))**um)**2.D0
              temp=1.D0/temp
              kr_s=0.5d0*kr(1)/se+2.d0*sqrt(se)*(1.d0-(1.d0-temp)**um)* &
                     (1.d0-temp)**(um-1.d0)*temp/se
              krg_s=-0.5D0/(1.D0-se)*kr(2) -2.D0*sqrt(1.D0-se)*((1.D0-se**(1.D0/um))**um)&
               *((1.D0-se**(1.D0/um))**(um-1.D0)) * (se**(1.D0/um-1.D0))
        ser=(sw-swir)/(1.D0-swir)
              upc=upc+(ser-se)*upc_s
              kr(1)=kr(1)+(ser-se)*kr_s
              kr(2)=kr(2)+ (ser-se)*krg_s
              
        !kr(2)=1.D0-kr(1)
            else
              upc=upc-upc_s*se
              kr(1)=0.D0
              kr(2)=1.D0
            end if
  !         print *,'in pckr nond ',sw,se,upc,kr
          end if
 !        print *,'in pckr  ',um,se,(temp-1.D0)
         if(upc > pcmax) upc=pcmax  

       case(2) !Brooks-Corey
       
          lam=pckr_lambda
          ala=pckr_alpha
       !  swir=pckr_swir

          if(sw>(1.05D0*swir))then
            se=(sw-swir)/(sw0-swir)
            upc=se**(-1.D0/lam)/ala
            kr(1)=se**(2.d0/lam+3.d0)
            !kr(2)=1.D0- kr(1)
            kr(2)=(1.D0 - se)**2.D0 * (1.D0 -se**(2.D0/lam +1.D0)) 
      else   ! use linear extropolation
            se=(0.05d0*swir)/(1.D0-swir)
            upc=se**(-1.D0/lam)/ala
            upc_s=-upc/se/lam
            if(sw>swir)then
              kr(1)=se**(2.D0/lam+3.d0)
              kr_s=(2.d0/lam+3.d0)*kr(1)/se
              krg_s = -2.D0*kr(2)/(1.D0-se) -(2.D0+lam)/lam*(1.D0-se)**2.D0*(se**(2.D0/lam)) 
        ser=(sw-swir)/(1.D0-swir)
              upc=upc+(ser-se)*upc_s
              kr(1)=kr(1)+(ser-se)*kr_s
        kr(2)=kr(2)+(ser-se)*krg_s
              !kr(2)=1.D0-kr(1)
            else
              upc=upc-upc_s*se
              kr(1)=0.D0
              kr(2)=1.D0
            end if
          end if
           
        case(3) !linear intropolation ,need pcmax, assign krmax=1.
        
           if(sw>swir)then
              se=(sw-swir)/(sw0-swir)
              upc=pcmax*(1.D0-se)
              kr(1)=se
              kr(2)=1.D0 - kr(1)
           else
              upc=pcmax
              kr(1)=0.d0
              kr(2)=1.d0
           end if


       case(4)  ! po model with gas phase residual
        
      if(sw>1.D0) sw=1.D0
        if(sw<0.D-0) sw= 1.D-5
      if(sw> pckr_sat_water_cut)then
            upc=0.D0; kr(1)=1.0d0; kr(2)=0.D0
      else
            ala=pckr_alpha
            betac=pckr_betac 
            um = pckr_m
            uum = 1.D0/um 
            un=1.D0/(1.D0-um)

            se=sw
            temp=se**uum
            upc=(1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
       
       if(upc>pcmax) upc=pcmax
            se=(sw-swir)/(1.D0-swir)
             st= 1.D0
            if(sw>=swir)then
        kr(1)= sqrt(se)*(1.D0-(1.D0-se**uum)**um)**2.D0
!             kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
              kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**pckr_pwr
            else
!        if(se<=0.D0) se= 1.D-7
        kr(1)=0.D0
!             kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
              kr(2)= sqrt(st-se)!*((1.D0-se**uum)**um)**pckr_pwr
               
            endif
          endif
        
          if(kr(1)<0.D0) kr(1)=0.D0!;kr(2)=1.D0
          if(kr(1)>1.D0) kr(1)=1.D0!;kr(2)=0.D0
          if(kr(2)<0.D0) kr(2)=0.D0!;kr(1)=1.D0
          if(kr(2)>1.D0) kr(2)=1.D0!;kr(1)=0.D0
        
        end select

        pc(1)=upc;  pc(2)=0.d0
        
!       if(sw<pckr_sat_water_cut) print *, sg,pc,kr
!      print *,'Ven ::',pc,kr

      return

    end subroutine pflow_pckr_noderiv_org

!_________________________________________________________________


subroutine pflow_pckr_richard(ipckrtype,pckr_swir,pckr_lambda, &
                 pckr_alpha,pckr_m,pckr_pcmax,sw,pc,kr,pckr_beta,pckr_pwr) 
      
      implicit none 


      integer ipckrtype
   !formation type, in pflow should be refered by grid%icap_loc
      real*8 :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax
      real*8 :: pckr_beta,pckr_pwr
      real*8 :  sw
      real*8 :: pc(*),kr(*)
       
      real*8 :: se,swir,sw0,lam,ala,um,un,upc,upc_s,kr_s, krg_s
      real*8 :: temp,ser,pcmax
      real*8 :: uum,pckr_betac,betac,st
     
     
    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
           
      sw0=1.d0
      pcmax=pckr_pcmax
      swir=pckr_swir
      select case(ipckrtype)

      case(1) ! van Gennuchten
          ala=pckr_alpha
      !    swir=pckr_swir
          um=pckr_m
          un=1.D0/(1.D0-um)
          se=(1.D0 + ala* pc(1)**un)**(-um)
          kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.d0
          sw = se* (sw0-swir) + swir
    
        case(2) !Brooks-Corey
      
          lam=pckr_lambda
          ala=pckr_alpha
       !  swir=pckr_swir

          se = (ala * pc(1))**(-lam)
          kr(1)=se**(2.d0/lam+3.d0)
          sw = se* (sw0-swir) + swir
        case(3) !linear intropolation ,need pcmax, assign krmax=1.
           se =1.D0 - pc(1)/ pcmax
           kr(1)= se
           sw = se* (sw0-swir) + swir
       end select pflow_pckr_richard
      
end subroutine 

!------------------------------------------------------------------------

    subroutine pflow_pckr(ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
              pckr_m ,pckr_pcmax,sg,pc,pc_s,kr,kr_s,pckr_beta,pckr_pwr) 
       
     
      integer :: ipckrtype
      real*8 :: sg
      real*8 :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_pwr
      real*8 :: pc(1:2),pc_s(1:2),kr(1:2),kr_s(1:2)
      real*8 :: pckr_beta
     
      real*8 :: sw,se,swir,sw0,lam,ala,um,un,upc,upc_s
      real*8 :: temp,pcmax,ser
      real*8 :: uum,pckr_betac,betac,st

    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
      sw=1.D0-sg
      swir=pckr_swir
      sw0=1.d0
      pcmax=pckr_pcmax
       
!     print *,'pflow_pckr: ',ipckrtype,sg,pckr_pcmax

      select case(ipckrtype)

      case(1) ! van Gennuchten
        ala=pckr_alpha
        um=pckr_m
        un=1.D0/(1.D0-um)

        if(sw>pckr_sat_water_cut)then
            upc=0.D0; upc_s=0.D0
            kr(1)=1.d0;  kr(2)=0.d0;
            kr_s(1)=0.d0; kr_s(2)=0.d0;
        elseif(sw>(1.05D0*swir))then
            se=(sw-swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
            temp=1.D0/temp
            kr(1)=sqrt(se)*(1.D0-(1.D0-temp)**um)**2.D0
            kr_s(1)=0.5d0*kr(1)/se+2.D0*sqrt(se)*(1.d0-(1.d0-temp)**um)* &
                  (1.d0-temp)**(um-1.d0)*temp/se
            kr(2)=1.D0-kr(1)
            kr_s(2)= -kr_s(1)
!           print *,'in pckr  ',um, sw , ala, se,upc,kr 
        else  ! use linear extropolation
            se=(0.05d0*swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
            if(sw>swir)then
              temp=1.D0/temp
              kr(1)=sqrt(se)*(1.D0-(1.D0-temp)**um)**2.D0
              kr_s(1)=0.5d0*kr(1)/se+2.d0*sqrt(se)*(1.D0-(1.D0-temp)**um)* &
                     (1.D0-temp)**(um-1.d0)*temp/se
              ser=(sw-swir)/(1.D0-swir)
              upc=upc+(ser-se)*upc_s
!             print *,se,ser,kr(1),kr_s(1)
              kr(1)=kr(1)+(ser-se)*kr_s(1)
              kr(2)=1.D0-kr(1)
              kr_s(2)= -kr_s(1)
            else
              upc=upc-upc_s*se
              upc_s=0.D0
              kr(1)=0.d0
              kr_s(1)=0.d0
              kr(2)=1.d0
              kr_s(2)=0.d0
            end if
          end if

        pc(1)=upc;  pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp


       case(2) !Brooks-Corey
          lam=pckr_lambda
          ala=pckr_alpha

          if(sw>(1.05D0*swir))then
             se=(sw-swir)/(sw0-swir)
             upc=se**(-1.D0/lam)/ala
             upc_s=-upc/se/lam
             kr(1)=se**(2.d0/lam+3.d0)
             kr_s(1)=(2.d0/lam+3.d0)*kr(1)/se
             kr(2)=1.D0-kr(1)
             kr_s(1)= -kr_s(1)
           else   ! use linear extropolation
              se=(0.05d0*swir)/(1.D0-swir)
              upc=se**(-1.D0/lam)/ala
              upc_s=-upc/se/lam
              if(sw>swir)then
                 kr(1)=se**(2.d0/lam+3.d0)
                 kr_s(1)=(2.d0/lam+3.d0)*kr(1)/se
                 ser=(sw-swir)/(1.D0-swir)
                 upc=upc+(ser-se)*upc_s
                 kr(1)=kr(1)+(ser-se)*kr_s(1)
                 kr(2)=1.D0-kr(1)
                 kr_s(2)=-kr_s(1)
              else
                 upc=upc-upc_s*se
                 kr(1)=0.D0
                 kr_s(1)=0.D0
                 kr(2)=1.D0
                 kr_s(2)=0.D0
              end if
           end if
        pc(1)=upc;  pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp


        case(3) !linear intropolation ,need pcmax, assign krmax=1.
         
           if(sw>swir)then
              se=(sw-swir)/(sw0-swir)
              upc=pcmax*(1.D0-se)
              upc_s= - pcmax
              kr(1)=se
              kr(2)=1.D0 - kr(1)
              kr_s(1)=1.D0
              kr_s(2)=-1.D0
           else
              upc=pcmax
              upc_s=0.D0
              kr(1)=0.d0
              kr(2)=1.d0
              kr_s(1)=0.d0
              kr_s(2)=0.d0
           end if
       
        pc(1)=upc;  pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp
       
       
      case(4)  ! po model with gas phase residual
        
      if(sw>1.D0) sw=1.D0
        if(sw<0.D-0) sw= 0.D-0
      
      if(sw> pckr_sat_water_cut)then
             upc=0.D0;  kr(1)=1.0;  kr(2)=0.D0;
       upc_s=0.D0; kr_s(1)=0.D0; kr(2)=0.D0;
      else
        ala=pckr_alpha
        betac=pckr_betac 
        um = pckr_m
        uum = 1.D0/um 
        un=1.D0/(1.D0-um)

        se=sw  
        temp=se**uum
        upc=(1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
        upc_s= -((1.D0/temp -1.D0)**(-um))*(se**(-uum-1.D0))*uum*(1.D0-um)&
           /ala/betac
          
        se=(sw-swir)/(1.D0-swir)
        st= 1.D0
              if(sw>=swir)then
          kr(1)= sqrt(se)*(1.D0-(1.D0-se**uum)**um)**2.D0
!               kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um )**7.D0
          kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um )**pckr_pwr
       
          kr_s(1) = 0.5D0 / sqrt(se)*(1.D0 - (1.D0 - se**uum)**um)**2.D0 &
            + 2.D0 * sqrt(se) * (1.D0- (1.D0-se**uum)**um) &
            * ((1.D0 - se**uum)**(um-1.D0)) * (Se**(uum-1.d0))
!               kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 -  se**uum)**(7.D0*um) &
!                     - 7.D0 * sqrt(1.D0- se) * (1.D0 - se**uum)**(7.D0*um-1.D0) &
!                     * se**(uum -1.D0)    
                kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0-se**uum)**(pckr_pwr*um) &
                - pckr_pwr * sqrt(1.D0-se) * (1.D0-se**uum)**(pckr_pwr*um-1.D0) &
                * se**(uum-1.D0)    
              else  
          kr(1)=0.D0
                kr_s(1)=0.D0

!               kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um)**7.D0
                kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um)**pckr_pwr
!               kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 - se**uum)**(7.D0*um) &
!               - 7.D0 * sqrt(1.D0- se) * (1.D0 - se**uum)**(7.D0*um-1.D0) &
!               * se**(uum -1.D0)    
                kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 - se**uum)**(pckr_pwr*um) &
                - pckr_pwr * sqrt(1.D0- se) * (1.D0 - se**uum)**(pckr_pwr*um-1.D0) &
                * se**(uum -1.D0)    
        endif
            endif      
    
            if(kr(1)<0.D0) kr(1)=0.D0; kr_s(1)=0.D0!; kr(2)=1.D0; kr_s(2)=0.D0
            if(kr(1)>1.D0) kr(1)=1.D0; kr_s(1)=0.D0!;kr(2)=0.D0; kr_s(2)=0.D0
            if(kr(2)<0.D0) kr(2)=0.D0; kr_s(2)=0.D0!;kr(1)=1.D0; kr_s(1)=0.D0
            if(kr(2)>1.D0) kr(2)=1.D0; kr_s(2)=0.D0!;kr(1)=0.D0; kr_s(1)=0.D0

            
        pc(1)=upc;  pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) 

    !if(pc(1)>pcmax) print *, 'pckr4: ',sg,pc,kr,pc_s,kr_s
        !if(sw<pckr_sat_water_cut) print *, sg,pc,kr,pc_s,kr_s
       
        end select

      return

    end subroutine pflow_pckr


! -----------------------
! New pckr routine for mph, vadose, owg

      subroutine pflow_pckr_noderiv(ipckrtype,pckr_sir,pckr_lambda, &
                 pckr_alpha,pckr_m,pckr_pcmax,sg,pc,kr,pckr_beta,pckr_pwr) 
      
      implicit none 


      integer ipckrtype
   !formation type, in pflow should be refered by grid%icap_loc
      real*8 :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax
      real*8 :: pckr_beta,pckr_pwr
      real*8 :: sg
      real*8 :: pc(1:2),kr(1:2)
       
      real*8 :: se,swir,sw0,lam,ala,um,un,upc,upc_s,kr_s, krg_s
      real*8 :: temp,ser,pcmax,sw
      real*8 :: uum,pckr_betac,betac,st
     
     
    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
      sw=1.D0-sg
     
      sw0=1.d0
      pcmax=pckr_pcmax
      swir=pckr_sir(1)
      select case(ipckrtype)

      case(1) ! van Gennuchten
          ala=pckr_alpha
      !    swir=pckr_swir
          um=pckr_m
          un=1.D0/(1.D0-um)
          if(sw> pckr_sat_water_cut)then
            upc=0.D0; kr(1)=1.d0; kr(2)=0.d0;
          elseif(sw>(1.05D0*swir))then
            se=(sw-swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.d0
            !kr(2)=1.D0-kr(1)
      kr(2)=sqrt(1.D0 - se)*((1.D0-se**(1.D0/um))**um)**2.D0
       !    print *,'in pckr nond ',sw,se,upc,kr
          else  ! use linear extropolation
            se=(0.05D0*swir)/(1.D0-swir)
            temp=se**(-1.D0/um)
            upc=(temp-1.D0)**(1.d0/un)/ala
            upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
            if(sw>swir)then
              kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.D0
        kr(2)=sqrt(1.D0 - se)*((1.D0-se**(1.D0/um))**um)**2.D0
              temp=1.D0/temp
              kr_s=0.5d0*kr(1)/se+2.d0*sqrt(se)*(1.d0-(1.d0-temp)**um)* &
                     (1.d0-temp)**(um-1.d0)*temp/se
              krg_s=-0.5D0/(1.D0-se)*kr(2) -2.D0*sqrt(1.D0-se)*((1.D0-se**(1.D0/um))**um)&
               *((1.D0-se**(1.D0/um))**(um-1.D0)) * (se**(1.D0/um-1.D0))
        ser=(sw-swir)/(1.D0-swir)
              upc=upc+(ser-se)*upc_s
              kr(1)=kr(1)+(ser-se)*kr_s
              kr(2)=kr(2)+ (ser-se)*krg_s
              
        !kr(2)=1.D0-kr(1)
            else
              upc=upc-upc_s*se
              kr(1)=0.D0
              kr(2)=1.D0
            end if
  !         print *,'in pckr nond ',sw,se,upc,kr
          end if
 !        print *,'in pckr  ',um,se,(temp-1.D0)
         if(upc > pcmax) upc=pcmax  

       case(2) !Brooks-Corey
       
          lam=pckr_lambda
          ala=pckr_alpha
       !  swir=pckr_swir

          if(sw>(1.05D0*swir))then
            se=(sw-swir)/(sw0-swir)
            upc=se**(-1.D0/lam)/ala
            kr(1)=se**(2.d0/lam+3.d0)
            !kr(2)=1.D0- kr(1)
            kr(2)=(1.D0 - se)**2.D0 * (1.D0 -se**(2.D0/lam +1.D0)) 
      else   ! use linear extropolation
            se=(0.05d0*swir)/(1.D0-swir)
            upc=se**(-1.D0/lam)/ala
            upc_s=-upc/se/lam
            if(sw>swir)then
              kr(1)=se**(2.D0/lam+3.d0)
              kr_s=(2.d0/lam+3.d0)*kr(1)/se
              krg_s = -2.D0*kr(2)/(1.D0-se) -(2.D0+lam)/lam*(1.D0-se)**2.D0*(se**(2.D0/lam)) 
        ser=(sw-swir)/(1.D0-swir)
              upc=upc+(ser-se)*upc_s
              kr(1)=kr(1)+(ser-se)*kr_s
        kr(2)=kr(2)+(ser-se)*krg_s
              !kr(2)=1.D0-kr(1)
            else
              upc=upc-upc_s*se
              kr(1)=0.D0
              kr(2)=1.D0
            end if
          end if
           
        case(3) !linear intropolation ,need pcmax, assign krmax=1.
        
           if(sw>swir)then
              se=(sw-swir)/(sw0-swir)
              upc=pcmax*(1.D0-se)
              kr(1)=se
              kr(2)=1.D0 - kr(1)
           else
              upc=pcmax
              kr(1)=0.d0
              kr(2)=1.d0
           end if


       case(4)  ! po model with gas phase residual
        
      if(sw>1.D0) sw=1.D0
        if(sw<0.D-0) sw= 1.D-5
      if(sw> pckr_sat_water_cut)then
            upc=0.D0; kr(1)=1.0d0; kr(2)=0.D0
      else
            ala=pckr_alpha
            betac=pckr_betac 
            um = pckr_m
            uum = 1.D0/um 
            un=1.D0/(1.D0-um)

            se=sw
            temp=se**uum
            upc=(1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
       
       if(upc>pcmax) upc=pcmax
            se=(sw-swir)/(1.D0-swir)
             st= 1.D0
            if(sw>=swir)then
        kr(1)= sqrt(se)*(1.D0-(1.D0-se**uum)**um)**2.D0
!             kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
              kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**pckr_pwr
            else
!        if(se<=0.D0) se= 1.D-7
        kr(1)=0.D0
!             kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
              kr(2)= sqrt(st-se)!*((1.D0-se**uum)**um)**pckr_pwr
               
            endif
          endif
        
          if(kr(1)<0.D0) kr(1)=0.D0!;kr(2)=1.D0
          if(kr(1)>1.D0) kr(1)=1.D0!;kr(2)=0.D0
          if(kr(2)<0.D0) kr(2)=0.D0!;kr(1)=1.D0
          if(kr(2)>1.D0) kr(2)=1.D0!;kr(1)=0.D0
        
        end select

        pc(1)=upc;  pc(2)=0.d0
        
!       if(sw<pckr_sat_water_cut) print *, sg,pc,kr
!      print *,'Ven ::',pc,kr

      return

    end subroutine pflow_pckr_noderiv

!_________________________________________________________________





  end module pckr_module
