

   subroutine oil_pckr_noderiv(ipckrtype,pckr_swir, pckr_soir, pcpckr_lambda, &
                 pckr_alpha,pckr_m,pckr_pcmax,s_w,s_o,pc,kr,pckr_beta,pckr_pwr) 
      
    integer :: ipckrtype
      real*8 :: s_w,s_o
      real*8 :: pckr_swir,pckr_soir, pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_pwr
      real*8 :: pc(1:3),kr(1:3)
      real*8 :: pckr_beta
	   
      real*8 :: sw,se,swir,sw0,lam,ala,um,un,upc,upc_s
      real*8 :: temp,pcmax,ser
      real*8 :: uum,pckr_betac,betac,st


     real*8, save :: alaph_ow= 1.50D-3, alaph_ao=2.64D-4
    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
      sw=1.D0-sg
      swir=pckr_swir
      sw0=1.d0
      pcmax=pckr_pcmax

	   
    pc(2)=0.D0
			 		     
	select case(ipckrtype)	     
	   case(11) ! Parker 1987
	      m_r = 1.D0/pckr_m
	      pckr_un=1.D0/(1.D0- pckr_m)
		  
		  
		  se_w = (s_w - pckr_swir)/( 1.D0 - pckr_swir)
		  se_t = (s_w + s_o - pckr_swir)/( 1.D0 - pckr_swir) 
		  if(se_w<1D-5)then
		    pc(1)=pcmax
		   else
			 pc(1)=((se_w)**pckr_m-1D0)**(1.D0/pckr_un)/alaph_ow
		     if(pc(1)>pcmax) pc(1)=pcmax
		   endif
		   if(se_t<1D-5)then
		    pc(3)=pcmax
		   else
			 pc(3)=((se_t)**pckr_m-1D0)**(1.D0/pckr_un)/alaph_ao
		     if(pc(3)>pcmax) pc(1)=pcmax
		   endif

		   
		  if (se_w<1D-5)then ! no water phase
		     kr(1)=0.D0
			  if( se_t<1D-5)then !no oil phase
			      kr(2)=1.D0; kr(3)=0.D0
			  elseif((1.D0-se_t)<1D-5)then !no gas phase
			    kr(2)=0.D0; kr(3)=1.D0
			  else         ! oil and gas phase
			    kr(2)= (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)	 
			    kr(3)=  (se_t-se_w)**5D-1 *(1.D0  -(1.D0 - se_t**m_r)**pckr_m)**2.D0
		      endif
		  else    ! have water
		     kr(1) = se_w **.5D0 * (1.D0 - (1.D0 - se_w**m_r)**pckr_m)**2.D0
			 if( (se_t-se_w)<1D-5)then ! no oil phase
			    kr(2)= (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)
                kr(3)=0.D0
			 elseif((1.D0-se_t)<1D-5)then !no gas phase
			    kr(2)=0.D0
				if((se_t-se_w)<1D-5)then ! no oil phase 
				  kr(3)=0.D0
				else
				  kr(3)= (se_t-se_w)**5D-1 *((1.D0 - se_w**m_r)**pckr_m -(1.D0 - se_t**m_r)**pckr_m)**2.D0
                endif
			 else	 									   				  		    		  
		  	   kr(1) = se_w **.5D0 * (1.D0 - (1.D0 - se_w**m_r)**pckr_m)**2.D0
			   kr(2) = (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)    
	           kr(3) = (se_t-se_w)**5D-1 *((1.D0 - se_w**m_r)**pckr_m -(1.D0 - se_t**m_r)**pckr_m)**2.D0
	         endif 
		  endif 
	     print *, s_w, s_o, se_w, se_t, pc,kr
	     case(12) ! Stone I 
	   
	   
	   case(13) ! Stone II		   	  
	  end select
	end subroutine oil_pckr_noderiv   


program main()

!use oil_pckr_module
implicit double precision(a-h,o-z)

real*8 pc(1:3), kr(1:3)
     k=11
     pckr_m=0.75D0
     ambda= pckr_m/ (1.D0- pckr_m)
       
     pckr_alaph=1D-3  
     pckr_pcmax=1D8
     
     pckr_beta=0.D0
     pckr_pwr=1.D0

  open(60, file="pckr.dat", action="write", status="new")
 
 do i=10,95
   do j=0,100-i
     sw=real(i)/100D0
     so=real(j)/100D0
     
      
     call oil_pckr_noderiv(k , 0.1D0, 0.2D0,ambda,pckr_alpha,pckr_m,pckr_pcmax,sw,so,pc,kr,pckr_beta,pckr_pwr)
                 
     write(60,'(8g12.4)') sw,so,pc,kr
     enddo
     enddo
    close(60)
    end program main