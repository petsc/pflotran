 module IMS_patch_module
  use pflow_gridtype_module

  private
! Cutoff parameters
  real*8, parameter :: formeps   = 0.D-5
  real*8, parameter :: eps       = 0.D-5
  real*8, parameter :: floweps   = 0.D-20
  real*8, parameter :: satcuteps = 0.D-5
  real*8, parameter :: dfac = 1.D-8
  real*8, pointer, save :: Resold_AR(:,:), Resold_FL(:,:)
 
! Contributions to residual from accumlation/source/Reaction, flux(include diffusion)
  
  public IMSRes_ARCont, IMSResidual_patch, IMSJacobin_patch

  contains 


  subroutine IMSRes_ARCont(node_no, var_node,por,vol,rock_dencpr, locpat,grid, Res_AR,ireac,ierr)
  implicit none
  integer node_no
  integer, optional:: ireac,ierr
  type(pflowGrid), intent(in) :: grid
  type(pflow_localpatch_info) :: locpat
  real*8, target:: var_node(1:grid%size_var_use)
  real*8 Res_AR(1:grid%ndof) 
  real*8 vol,por,rock_dencpr
	   
  real*8, pointer :: temp, pre_ref   ! 1 dof
  real*8, pointer :: sat(:), density(:), amw(:), h(:), u(:), pc(:), kvr(:)         ! nphase dof
    
  integer :: ibase, m,np, iireac=1
  real*8 pvol,mol(grid%nphase),eng
  
  if(present(ireac)) iireac=ireac
  pvol=vol*por
  
  ibase=1;                 temp=>var_node(ibase)
  ibase=ibase+1;           pre_ref=>var_node(ibase)
  ibase=ibase+1;           sat=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr=>var_node(ibase:ibase+grid%nphase-1)

  !sumation of component
  mol=0.D0; eng=0.D0
  do np = 1, grid%nphase
   if(sat(np)>1.D-11)then  
   
     mol(np) = mol(np) + sat(np)*density(np)
     !eng = eng + density(np)*u(np) *sat(np)
  endif
  enddo

  mol = mol * pvol
  !eng = eng * pvol + (1.D0 - por)* vol * rock_dencpr * temp 
  
   Res_AR(1:grid%nphase)=mol(:)
   !Res_AR(grid%ndof)=eng
  nullify(temp, pre_ref, sat, density, pc, kvr) 	    
  end subroutine  IMSRes_ARCont


 subroutine IMSRes_FLCont(nconn_no,area, &
         var_node1,por1,tor1,sir1,dd1,perm1,Dk1,&
		 var_node2,por2,tor2,sir2,dd2,perm2,Dk2,&
		 locpat, grid, vv_darcy,Res_FL)
 implicit none
  integer nconn_no
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info) :: locpat
  real*8 sir1(1:grid%nphase),sir2(1:grid%nphase)
  real*8, target:: var_node1(1:2+4*grid%nphase)
  real*8, target:: var_node2(1:2+4*grid%nphase)
  real*8 por1,por2,tor1,tor2,perm1,perm2,Dk1,Dk2,dd1,dd2
  real*8 vv_darcy(grid%nphase),area
  real*8 Res_FL(1:grid%ndof) 
	   
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  
  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
     
  
  integer ibase, m,np, ind
  real*8  fluxm(grid%nphase),fluxe, v_darcy,q
  real*8  ukvr, DK,Dq, diffdp
  real*8 upweight,density_ave,cond, gravity, dphi
  

  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
						   
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
						   
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+grid%nphase-1)
						   sat2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density1=>var_node1(ibase:ibase+grid%nphase-1)
                           density2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc1=>var_node1(ibase:ibase+grid%nphase-1)
                           pc2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr1=>var_node1(ibase:ibase+grid%nphase-1)
                           kvr2=>var_node2(ibase:ibase+grid%nphase-1)
  
  !print *,' FLcont got pointers' ,var_node1,var_node2,sir1,sir2
  !print *,' tmp=',temp1,temp2
  !print *,'diff=',diff1,diff2
   
  Dq=(perm1 * perm2)/(dd1*perm2 + dd2*perm1)
  diffdp= (por1 *tor1 * por2*tor2) / (dd2*por1*tor1 + dd1*por2*tor2)*area
  
  fluxm=0.D0
  fluxe=0.D0
  vv_darcy=0.D0  
  
  do np=1, grid%nphase

! Flow term
    if ((sat1(np) > sir1(np)) .or. (sat2(np) > sir2(np)))then
	  
	  upweight=.5D0
	  if(sat1(np) <eps) then 
         upweight=0.d0
	    else if(sat2(np) <eps) then 
         upweight=1.d0
      endif
	  density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
	  
	  gravity = (upweight*density1(np) + (1.D0-upweight)*density2(np)) &
              * grid%gravity * locpat%delz(nconn_no)
			  
	  dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
!	  print *,'FLcont  dp',dphi
	! note uxmol only contains one phase xmol
	  if(dphi>=0.D0)then
	     ukvr=kvr1(np)
		else
		 ukvr=kvr2(np)
	  endif 	   
     
	 ! print *,'FLcont  uxmol',uxmol
	  if(ukvr>floweps)then
         v_darcy= Dq * ukvr * dphi
         !grid%vvl_loc(nconn_no) = v_darcy
		 vv_darcy(np)=v_darcy
		 
	     q=v_darcy * area

         fluxm(np)=fluxm(np) + q*density_ave																														
																																																																																								    
          ! fluxe = fluxe + q*density_ave*uh 
       endif
   endif 
 !  print *,' FLcont end flow',np
 enddo
     
! conduction term
        
        !Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)
		!cond=Dk*area*(temp1-temp2) 
        !fluxe=fluxe + cond
   !      print *,' FLcont heat cond', Dk, cond
 Res_FL(1:grid%nphase)=fluxm(:) * grid%dt
 !Res_FL(grid%ndof)=fluxe * grid%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 !												   2 R = R - Res_FL	
 !print *,'end FLcont'
 
  nullify(temp1, pre_ref1, sat1, density1, pc1,kvr1) 	    
  nullify(temp2, pre_ref2, sat2, density2, pc2,kvr2) 	    
 end subroutine IMSRes_FLCont

 subroutine IMSRes_FLBCCont(nbc_no,area, &
             var_node1,var_node2,por2,tor2,sir2,dd1,perm2,Dk2,&
			 locpat, grid, vv_darcy,Res_FL)
 ! Notice : index 1 stands for BC node
   implicit none
  
   integer nbc_no
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info) :: locpat
  real*8 dd1, sir2(1:grid%nphase)
  real*8, target:: var_node1(1:2+4*grid%nphase)
  real*8, target:: var_node2(1:2+4*grid%nphase)
  real*8 por2,perm2,Dk2,tor2
  real*8 vv_darcy(grid%nphase), area
  real*8 Res_FL(1:grid%ndof) 
	   
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), pc1(:), kvr1(:)         ! nphase dof

  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), pc2(:), kvr2(:)         ! nphase dof

  
  integer ibase, m,np, ind, ibc,j
  real*8  fluxm(grid%nphase),fluxe, v_darcy,q
  real*8 uh, ukvr,diff,diffdp, DK,Dq
  real*8 upweight,density_ave,cond,gravity, dphi

  
  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+grid%nphase-1)
						   sat2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density1=>var_node1(ibase:ibase+grid%nphase-1)
                           density2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc1=>var_node1(ibase:ibase+grid%nphase-1)
                           pc2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr1=>var_node1(ibase:ibase+grid%nphase-1)
                           kvr2=>var_node2(ibase:ibase+grid%nphase-1)


   ibc = locpat%ibconn(nbc_no)
   fluxm=0.D0; fluxe=0.D0
   vv_darcy=0.D0 
   
   select case (grid%ibndtyp(ibc))
     case(1) 
	   Dq= perm2 / dd1
        diffdp = por2*tor2/dd1*area
        ! Flow term
		do np =1,grid%nphase
         if ((sat1(np) > sir2(np)) .or. (sat2(np) > sir2(np)))then
	  
	     upweight=.5D0
	     if(sat1(np) <eps) then 
             upweight=0.d0
	        else if(sat2(np) <eps) then 
             upweight=1.d0
         endif
	  density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
	  
	  gravity = (upweight*density1(np)+ (1.D0-upweight)*density2(np)) &
              * grid%gravity * locpat%delzbc(nbc_no)
			  
	  dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
	  
	 
	  if(dphi>=0.D0)then
	      ukvr=kvr1(np)
		else
		  ukvr=kvr2(np)
	  endif 	   
     
	  if(ukvr*Dq>floweps)then
         v_darcy= Dq * ukvr * dphi
         !grid%vvl_loc(nbc_no) = v_darcy
		 vv_darcy(np)=v_darcy
		 
	     q=v_darcy * area
		 		 
		 fluxm(np)=fluxm(np) + q*density_ave																														
	     																																																																																									    
          ! fluxe = fluxe + q*density_ave*uh 
       endif
      endif 
   enddo
! conduction term
        
        !Dk =  Dk2 / dd1
		!cond=Dk*area*(temp1-temp2) 
        !fluxe=fluxe + cond
 
     Res_FL(1:grid%nphase)=fluxm(:)* grid%dt
    ! Res_FL(grid%ndof)=fluxe * grid%dt

  case(2)
    if((dabs(locpat%velocitybc(1,nbc_no))+dabs(locpat%velocitybc(2,nbc_no)))>floweps)then
      ! print *, 'FlowBC :', nbc_no,locpat%velocitybc(1,nbc_no),locpat%velocitybc(2,nbc_no)
	   do j=1,grid%nphase
          fluxm=0.D0; fluxe=0.D0
          v_darcy = locpat%velocitybc(j,nbc_no)
		  vv_darcy(j) = locpat%velocitybc(j,nbc_no)
!		  grid%vvbc(j+(nc-2)*grid%nphase)= grid%velocitybc(j,nc)
	    ! note different from 2 phase version

         if(v_darcy >0.d0)then 
             q = v_darcy * density1(j) * area
             !q = 0.d0
             !flux = flux - q
          !   fluxe = fluxe - q  * h1(j) 

                fluxm(j)=fluxm(j) + q 
  
           else 
              q =  v_darcy * density2(j) * area   
           !   fluxe = fluxe - q  * h2(j) 
			  fluxm(j)=fluxm(j) + q 
             
          endif 
   
	enddo
    endif
	Res_FL(1:grid%nphase)=fluxm(:)* grid%dt
   ! Res_FL(grid%ndof)=fluxe * grid%dt

   case(3)
     Dq = perm2/dd1 

     do np =1,grid%nphase
      if ((sat2(np) > sir2(np)))then
	  
	     density_ave = density2(np)  
	  
	     gravity = density2(np)* grid%gravity * locpat%delzbc(nbc_no)
	 		  
	     dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
	  
	     ukvr=kvr2(np)
	   
	   if(ukvr*Dq>floweps)then
         v_darcy= Dq * ukvr * dphi
    !     grid%vvbc(,nbc_no) = v_darcy
		 vv_darcy(np)=v_darcy
		 
	     q=v_darcy * area
		 		 
		
		   fluxm(np)=fluxm(np) + q*density_ave																														
	    																																																																																										    
          ! fluxe = fluxe + q*density_ave*uh 
        endif
     endif 
 
     enddo

    Res_FL(1:grid%nphase)=fluxm(:)* grid%dt
    !Res_FL(grid%ndof)=fluxe * grid%dt

  end select
   nullify(temp1, pre_ref1, sat1, density1, pc1,kvr1) 	    
  nullify(temp2, pre_ref2, sat2, density2, pc2,kvr2) 	    
 end  subroutine IMSRes_FLBCCont 



 subroutine IMSResidual_patch(xx,r_p,accum_p, locpat,grid,ierr)
    use translator_IMS_module
   
    implicit none
 

! Here the xx works as an array, contains information as xx_loc_p

    real*8, pointer :: xx(:)

! r still using local index, not ghosted 	
    real*8, intent(out) :: r_p(:)
    real*8, pointer ::accum_p(:)
! var is no longer a petsc vector    
	
   type(pflowGrid), intent(inout) :: grid
   type(pflow_localpatch_info) :: locpat
   
 
  integer :: ierr
  integer :: n, ng, nc, nr
  integer :: i, i1, i2, j, jn, jng, jm1, jm2, jmu
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2, p1, p2, t1, t2, c1, c2,&
             s1, s2
  integer :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
  integer :: i1_hencoeff, i2_hencoeff
  integer :: ibc  ! Index that specifies a boundary condition block
  integer, save ::icall
! real*8 :: term1, term2, term3

!  real*8, pointer :: Resold_AR(:,:), Resold_FL(:,:) 

  real*8, pointer :: yy_p(:),&
                 porosity_loc_p(:), volume_p(:), &
                 phis_p(:), tor_loc_p(:),&
                 perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
                 vl_p(:), var_p(:),var_loc_p(:) 
                          
               
  real*8, pointer :: pc_p(:), pc_loc_p(:),kvr_p(:), kvr_loc_p(:)

  real*8, pointer :: icap_p(:), icap_loc_p(:), ithrm_loc_p(:),ithrm_p(:)
  real*8, pointer :: deltax(:)



  integer :: iicap, index_var_begin, index_var_end,iicap1,iicap2,np
  integer :: ichange
  real*8 :: dd1, dd2, eng, cond, den, temp,&
            eengl,eengg, &
            fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl,&
            fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q,  &
            v_darcy,hflx,pvoldt, voldt, accum, pvol
  real*8 :: dd, f1, f2, ff, por1, por2, perm1, perm2
  real*8 :: Dphi,D0
  real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol,density_ave,dif(grid%nphase)
  real*8 :: tsrc1, qsrc1(grid%nphase), maxqsrc
  real*8 :: cw,cw1,cw2, xxlw,xxla,xxgw,xxga
  real*8 :: upweight
  real*8 :: ukvr,uhh,uconc, tmp
  real*8 :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp, rho, visc
  real*8 :: Res(grid%ndof), vv_darcy(grid%nphase)
  data icall/0/
 
 
 
 temp=grid%tref


 ! call VecGetArrayF90(grid%var_loc,var_loc_p,ierr)
  var_loc_p => locpat%var
 ! call VecGetArrayF90(grid%yy,yy_p,ierr)
  yy_p => locpat%yy_p
 ! call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  porosity_loc_p => locpat%porosity_loc_p
 ! call VecGetArrayF90(grid%tor_loc, tor_loc_p, ierr)
  tor_loc_p => locpat%tor_loc_p
 ! call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  perm_xx_loc_p=> locpat%perm_xx_loc_p
 ! call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  perm_yy_loc_p=> locpat%perm_yy_loc_p
 ! call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  perm_zz_loc_p=> locpat%perm_zz_loc_p
 ! call VecGetArrayF90(grid%volume, volume_p, ierr)
  volume_p=> locpat%volume_p
 ! call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  ithrm_loc_p=> locpat%ithrm_loc_p
 ! call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
 icap_loc_p=> locpat%icap_loc_p
 ! call VecGetArrayF90(grid%vl, vl_p, ierr)
 vl_p=> locpat%vl_p
  !print *,' Res_patch: get pinter '

 !Resold_AR => locpat%Resold_AR
 !Resold_FL => locpat%Resold_FL
 
 if (icall==0)then
   
   allocate(Resold_AR(locpat%nlmax,grid%ndof))
   allocate(Resold_FL(locpat%nconn,grid%ndof))
   !print *,' Res_patch: get pinter '
   icall =1
 endif
 
   
!there is potential possiblity that the pertubate of p may change the direction of pflow.
! once that happens, code may crash, namely wrong matrix element 
 do ng = 1, locpat%ngmax 
  !ng=grid%nL2G(n)
  
    locpat%delx(1,ng)=xx((ng-1)*grid%ndof+1)*dfac*1e-3
  !  grid%delx(2,ng)=xx_loc_p((ng-1)*grid%ndof+2)*dfac
   
   do i=2, grid%ndof		
		  if(xx((ng-1)*grid%ndof+i) <=0.5)then
		   locpat%delx(i,ng) = dfac *xx((ng-1)*grid%ndof+i) 
		  else
		   locpat%delx(i,ng) = -dfac *xx((ng-1)*grid%ndof+i) 
		  endif 
		  
		  if(locpat%delx(i,ng) <1D-9 .and. locpat%delx(i,ng)>=0.D0)locpat%delx(i,ng) =1D-9
		  if(locpat%delx(i,ng) >-1D-9 .and. locpat%delx(i,ng)<0.D0)locpat%delx(i,ng) =-1D-9
		  
		  if((locpat%delx(i,ng)+xx((ng-1)*grid%ndof+i))>1.D0)then
			 locpat%delx(i,ng) = (1.D0-xx((ng-1)*grid%ndof+i))/1D5
		   endif
		  if((locpat%delx(i,ng)+xx((ng-1)*grid%ndof+i))<0.D0)then
			 locpat%delx(i,ng) = xx((ng-1)*grid%ndof+i)/1D5
		   endif
    enddo
  
  enddo
  
 

    !print *, 'IMS_Res: Got delx'
 ! call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)
!------------------------------------------------------ 




!-----  phase properities ---- last time step---
  do ng = 1, locpat%ngmax
    iicap = icap_loc_p(ng)

	call pri_var_trans_IMS_ninc(xx((ng-1)*grid%ndof+1:ng*grid%ndof),&
        temp, grid%scale,grid%nphase,&
        iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),&
		locpat%var((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),&
		ierr)



	if (grid%ideriv .eq. 1) then
      call pri_var_trans_ims_winc(xx((ng-1)*grid%ndof+1:ng*grid%ndof),&
	    locpat%delx(1:grid%ndof,ng), temp,&
        grid%scale,grid%nphase,&
        iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),&
        locpat%var((ng-1)*grid%size_var_node+grid%size_var_use+1:ng*grid%size_var_node),&
	    ierr)
    endif
	
!	print *,'var_p',n,iiphase, var_p((n-1)*grid%size_var_node+1:n*grid%size_var_node)					    
!   if(n < 5) print *,'pflow_2ph: ',n,grid%ideriv,grid%xxphi_co2(n)
  enddo

 

   
  

  ! notice:: here we assume porosity is constant
  
  ! for smallest change of code, I point pointers with old name style to the new one
 
 
  

  Resold_AR=0.D0; ResOld_FL=0.D0

!--------------------------------------------------------------------------
! Calculate accumulation term for interior and exterior nodes.
!--------------------------------------------------------------------------
  
  r_p = - accum_p
 ! in fortran, this statement will copy values, not duplicate pointer address 

  do n = 1, locpat%nlmax  ! For each local node do...
    ng = locpat%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
    index_var_begin=(ng-1)*grid%size_var_node+1
    index_var_end = index_var_begin -1 + grid%size_var_use
	  
    pvol = volume_p(n)*porosity_loc_p(ng)
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt
  
    i = ithrm_loc_p(ng)

    accum = 0.d0
    call IMSRes_ARCont(n, var_loc_p(index_var_begin: index_var_end),&
	  porosity_loc_p(ng),volume_p(n),grid%dencpr(i),locpat,grid, Res, 1,ierr)
   
    r_p(p1:p1+grid%ndof-1) = r_p(p1:p1+grid%ndof-1) + Res(1:grid%ndof)
    Resold_AR(n,1:grid%ndof)= Res(1:grid%ndof) 
  end do
  !print *,' Finished accum terms'

!************************************************************************
 ! add source/sink terms
 
  do nr = 1, grid%nblksrc
      
    kk1 = grid%k1src(nr) - locpat%nzs
    kk2 = grid%k2src(nr) - locpat%nzs
    jj1 = grid%j1src(nr) - locpat%nys
    jj2 = grid%j2src(nr) - locpat%nys
    ii1 = grid%i1src(nr) - locpat%nxs
    ii2 = grid%i2src(nr) - locpat%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(locpat%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(locpat%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(locpat%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1(:) = grid%qsrc(i,nr,:)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1(:) = f1*grid%qsrc(i,nr,:) + f2*grid%qsrc(i-1,nr,:)
        goto 10
      endif
    enddo
 10 continue
    
   !print *,'pflow2ph : ', grid%myrank,i,grid%timesrc(i,nr), &
   !grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
 	
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
	
	maxqsrc=0.D0
	do j=1, grid%nphase
	 if(dabs(qsrc1(j))>0.D0) maxqsrc =dabs(qsrc1(j))
	enddo 
	
    if (maxqsrc > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*locpat%nlx+(kk-1)*locpat%nlxy
            ng = locpat%nL2G(n)
            
		   
           !           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
            r_p((n-1)*grid%ndof +1:(n-1)*grid%ndof +grid%nphase) = &
			    r_p((n-1)*grid%ndof +1:(n-1)*grid%ndof +grid%nphase) - qsrc1(:) *grid%dt
          !  r_p(n*grid%ndof) = r_p(n*grid%ndof) - qsrc1*enth_src_h2o*grid%dt
            Resold_AR(n,:)= Resold_AR(n,:) - qsrc1(:)*grid%dt
		!	Resold_AR(n,grid%ndof)= Resold_AR(n,grid%ndof) - qsrc1 * enth_src_h2o*grid%dt

          enddo
        enddo
      enddo
    endif  
		
  enddo
  !print *,'finished source/sink term'
  

!*********************************************************************


 
! stop
!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------
 
  do nc = 1, locpat%nconn  ! For each interior connection...
    m1 = locpat%nd1(nc) ! ghosted
    m2 = locpat%nd2(nc)

    n1 = locpat%nG2L(m1) ! = zero for ghost nodes
    n2 = locpat%nG2L(m2) ! Ghost to local mapping   

    p1 = 1 + (n1-1)*grid%ndof 
    p2 = 1 + (n2-1)*grid%ndof
   
    dd1 = locpat%dist1(nc)
    dd2 = locpat%dist2(nc)
    
    ip1 = locpat%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = locpat%iperm2(nc)


    select case(ip1)
    case(1) 
       perm1 = perm_xx_loc_p(m1)
    case(2)
       perm1 = perm_yy_loc_p(m1)
    case(3)
       perm1 = perm_zz_loc_p(m1)
    end select
    
    select case(ip2)
    case(1) 
       perm2 = perm_xx_loc_p(m2)
    case(2)
       perm2 = perm_yy_loc_p(m2)
    case(3)
       perm2 = perm_zz_loc_p(m2)
    end select

       i1 = ithrm_loc_p(m1)
       i2 = ithrm_loc_p(m2)
	   iicap1=int(icap_loc_p(m1))
	   iicap2=int(icap_loc_p(m2))
	   
       D1 = grid%ckwet(i1)
       D2 = grid%ckwet(i2)

    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
    
    call IMSRes_FLCont(nc ,locpat%area(nc), &
         var_loc_p((m1-1)*grid%size_var_node+1:(m1-1)*grid%size_var_node+grid%size_var_use),&
		 porosity_loc_p(m1),tor_loc_p(m1),grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1,&
		 var_loc_p((m2-1)*grid%size_var_node+1:(m2-1)*grid%size_var_node+grid%size_var_use),&
		 porosity_loc_p(m2),tor_loc_p(m2),grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,&
		 locpat, grid, vv_darcy,Res)
     ! grid%vvl_loc(nc) = vv_darcy(1)
	 ! grid%vvg_loc(nc) = vv_darcy(2)  
    
	 
	!   use for print out of velocity **********************************************
	    if (n1 > 0) then               ! If the upstream node is not a ghost node...
        do np =1, grid%nphase 
		   vl_p(np+(ip1-1)*grid%nphase+3*grid%nphase*(n1-1)) = vv_darcy(np) 
!note: grid%vl is indexed according to n1 numbering 
        enddo
	 endif
    
	 
	  
	    
	 Resold_FL(nc,1:grid%ndof)= Res(1:grid%ndof) 
    
	if(n1>0)then
	   r_p(p1:p1+grid%ndof-1) = r_p(p1:p1+grid%ndof-1) + Res(1:grid%ndof)
     endif
   
    if(n2>0)then
	   r_p(p2:p2+grid%ndof-1) = r_p(p2:p2+grid%ndof-1) - Res(1:grid%ndof)
     endif


 end do
!  print *,'finished NC' 
 
!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', grid%myrank, grid%sgbc    
 
  do nc = 1, locpat%nconnbc

       m = locpat%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       ng = locpat%nL2G(m)

       if(ng<=0)then
         print *, "Wrong boundary node index... STOP!!!"
         stop
       end if
  
       p1 = 1 + (m-1) * grid%ndof
       
	   ibc = locpat%ibconn(nc)
       ip1 = locpat%ipermbc(nc)
	   

       i2 = ithrm_loc_p(ng)
       D2 = grid%ckwet(i2)


       select case(ip1)
         case(1)
           perm1 = perm_xx_loc_p(ng)
         case(2)
           perm1 = perm_yy_loc_p(ng)
         case(3)
           perm1 = perm_zz_loc_p(ng)
       end select

       select case(grid%ibndtyp(ibc))
          
          case(2)
          ! solve for pb from Darcy's law given qb /= 0
              locpat%xxbc(:,nc)= xx((ng-1)*grid%ndof+1: ng*grid%ndof)
                 case(3) 
     !         grid%xxbc((nc-1)*grid%ndof+1)=grid%pressurebc(2,ibc)
			  locpat%xxbc(2:grid%ndof,nc)= &
			           xx((ng-1)*grid%ndof+2: ng*grid%ndof)
			 ! grid%iphasebc(nc)=int(iphase_loc_p(ng))
		    end select

    
     iicap=int(icap_loc_p(ng))  
	
      call pri_var_trans_ims_ninc(locpat%xxbc(:,nc), temp,&
	     grid%scale,grid%nphase, &
      iicap,grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), & !use node's value
      grid%pcbetac(iicap),grid%pwrprm(iicap),&
      locpat%varbc(1:grid%size_var_use), ierr)
   
     
    call IMSRes_FLBCCont(nc,locpat%areabc(nc), &
             locpat%varbc(1:grid%size_var_use), &
			 var_loc_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),&
			 porosity_loc_p(ng),tor_loc_p(ng),grid%sir(1:grid%nphase,iicap),locpat%distbc(nc),&
			 perm1,D2, locpat,grid, vv_darcy,Res)
     locpat%vvlbc(nc) = vv_darcy(1)
	 locpat%vvgbc(nc) = vv_darcy(2) 
    r_p(p1:p1-1+grid%ndof)= r_p(p1:p1-1+grid%ndof) - Res(1:grid%ndof)
	ResOld_AR(m,1:grid%ndof) = ResOld_AR(m,1:grid%ndof) - Res(1:grid%ndof)
   
   
       !print *, ' boundary index', nc,ng,ibc,grid%ibndtyp(ibc)
       !print *,'        xxbc', locpat%iphasebc(nc), locpat%xxbc(:,nc),res
	   !print *, '       var', locpat%varbc
	!   print *, ' P  T   C   S  ', locpat%pressurebc(1,ibc),locpat%tempbc(ibc), &
    !                               locpat%concbc(ibc),locpat%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),grid%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,grid%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
!print *,grid%density_bc,grid%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc

 enddo
  !print *,'finished BC'
 


! print *,'finished IMSResidual'
! nullify(Resold_AR, Resold_FL)
 end subroutine IMSResidual_patch



subroutine IMSJacobin_patch(xx, Jac_elem, Jac_ind ,locpat,grid,ierr)
       
    use translator_IMS_module
   
    implicit none

    real*8, pointer :: xx(:)

    ! matrix block, 0 to 6 should be 3D non-zero entrys.!
	!               0 stands for diagnal element
	!               1 left.  2. right  3 up.  4 down.  5 front. 6 back.
	!               This is determined by the connection setup, does not really matter   
	
    type(pflowGrid), intent(inout) :: grid
	type(pflow_localpatch_info) :: locpat
   ! integer, intent(inout) :: flag
    real*8, intent(out) :: Jac_elem(1:locpat%nlmax,0:6,1:grid%ndof,1:grid%ndof)

	! non-zero matrix cloumn index, diagnal will be the same as row index
	! for petsc version, it should be global index 
	integer, intent(out) :: Jac_Ind(1:locpat%nlmax,0:6)

    integer :: ierr
    integer*4 :: n, ng, nc,nvar,neq,nr
    integer*4 :: i1, i2, j, jn, jng, jm1, jm2,jmu, i
    integer   :: kk,ii1,jj1,kk1,ii2,jj2,kk2  
    integer*4 :: m, m1, m2, mu, n1, n2, ip1, ip2 
    integer*4 :: p1,p2,t1,t2,c1,c2,s1,s2
    integer*4 :: ibc  ! Index that specifies a boundary condition block.
    real*8 ::  v_darcy, q

    real*8, pointer :: porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), phis_p(:),  tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
               
               
       


  real*8, pointer ::  icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  integer :: iicap,ii,jj,iiphas,iiphas1,iiphas2,iicap1,iicap2
  integer :: index_var_begin, index_var_end
  integer*4 ibc_hencoeff
  real*8 :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho,dddt,dddp,fg,dfgdp,&
            dfgdt,eng,dhdt,dhdp,visc,dvdt,dvdp
  real*8 :: cond, gravity,  acc,  vv_darcy(grid%nphase),&
            density_ave, voldt, pvoldt
  real*8 :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
            flux, fluxh, fluxv, ff
  real*8 :: tsrc1,qsrc1(grid%nphase)
  real*8 :: dd1, dd2, dd, f1, f2, den
! real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
  real*8 :: por1, por2, perm1, perm2
  real*8 :: qu_rate, p_vapor,sat_pressure_t
! real*8 :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
  real*8 :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu, sat_ave
  real*8 :: ra(1:grid%ndof,1:2*grid%ndof)  
  real*8 :: uhh, uconc, ukvr, tmp, temp
  real*8 :: upweight,m1weight,m2weight,mbweight,mnweight
  real*8 :: delxbc(1:grid%ndof)
  real*8 :: blkmat11(1:grid%ndof,1:grid%ndof), &
            blkmat12(1:grid%ndof,1:grid%ndof),&
			blkmat21(1:grid%ndof,1:grid%ndof),&
			blkmat22(1:grid%ndof,1:grid%ndof)
  real*8 :: ResInc(1:locpat%nlmax, 1:grid%ndof, 1:grid%ndof),res(1:grid%ndof)  
  real*8 :: max_dev, maxqsrc  
  integer  na1,na2,neighbor_ind1,neighbor_ind2,ndef
!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
 
 temp=grid%tref
Jac_elem =0.D0
Jac_ind = -1

 !print *,'*********** In Jacobian ********************** '

!  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
!   xx_loc_p => xx
!  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  porosity_loc_p=> locpat%porosity_loc_p
!  call VecGetArrayF90(grid%tor_loc, tor_loc_p, ierr)
  tor_loc_p=> locpat%tor_loc_p
!  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  perm_xx_loc_p=> locpat%perm_xx_loc_p
!  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  perm_yy_loc_p=> locpat%perm_yy_loc_p
!  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  perm_zz_loc_p=> locpat%perm_zz_loc_p
!  call VecGetArrayF90(grid%volume, volume_p, ierr)
  volume_p=> locpat%volume_p
!  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  ithrm_loc_p=> locpat%ithrm_loc_p
!  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  icap_loc_p=> locpat%icap_loc_p
!  call VecGetArrayF90(grid%var_loc, var_loc_p, ierr)
  var_loc_p=> locpat%var

 !Resold_AR => locpat%Resold_AR
 !Resold_FL => locpat%Resold_FL

 !print *,' In mph Jacobian ::  got pointers '
! ********************************************************************

! Accumulation terms

 ResInc=0.D0
 
 
  do n = 1, locpat%nlmax  ! For each local node do...
    ng = locpat%nL2G(n)   !get ghosted index
    
	voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

 !    iiphas=iphase_loc_p(ng)
 ! pressure equation    
   do nvar=1, grid%ndof
   
      index_var_begin=(ng-1)*grid%size_var_node+nvar*grid%size_var_use+1
      index_var_end = index_var_begin -1 + grid%size_var_use

       call IMSRes_ARCont(n, var_loc_p(index_var_begin : index_var_end),&
	      porosity_loc_p(ng),volume_p(n),grid%dencpr(int(ithrm_loc_p(ng))),&
	      locpat, grid, Res,1,ierr)
      
       ResInc(n,:,nvar) = ResInc(n,:,nvar) + Res(:)
   end do
 enddo
 !print *,' Mph Jaco Finished accum terms'
! Source / Sink term

   do nr = 1, grid%nblksrc
      
    kk1 = grid%k1src(nr) - locpat%nzs
    kk2 = grid%k2src(nr) - locpat%nzs
    jj1 = grid%j1src(nr) - locpat%nys
    jj2 = grid%j2src(nr) - locpat%nys
    ii1 = grid%i1src(nr) - locpat%nxs
    ii2 = grid%i2src(nr) - locpat%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(locpat%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(locpat%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(locpat%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1(:) = grid%qsrc(i,nr,:)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1(:) = f1*grid%qsrc(i,nr,:) + f2*grid%qsrc(i-1,nr,:)
        goto 10
      endif
    enddo
 10 continue
 
 
    maxqsrc=0.D0
	do j=1, grid%nphase
	 if(dabs(qsrc1(j))>0.D0) maxqsrc =dabs(qsrc1(j))
	enddo 

    
     if (maxqsrc > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*locpat%nlx+(kk-1)*locpat%nlxy
            ng = locpat%nL2G(n)
            
	        do nvar=1,grid%ndof 	   
            
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
            ResInc(n,1:grid%nphase,nvar)=  ResInc(n,1:grid%nphase,nvar) - qsrc1(:)*grid%dt
        !    ResInc(n,grid%ndof,nvar)=  ResInc(n,grid%ndof,nvar) - qsrc1*enth_src_h2o*grid%dt

			
			
			!           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*grid%fmwh2o, &
!           qsrc1
          enddo
		  enddo
        enddo
      enddo
    endif  
	
  enddo	
  
   !print *,' Mph Jaco Finished source terms'
! Contribution from BC
 do nc = 1, locpat%nconnbc

       m = locpat%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       ng = locpat%nL2G(m)

       if(ng<=0)then
         print *, "Wrong boundary node index... STOP!!!"
         stop
       end if
  
       p1 = 1 + (m-1) * grid%ndof
       
	   ibc = locpat%ibconn(nc)
       ip1 = locpat%ipermbc(nc)
	   
      
       i2 = ithrm_loc_p(ng)
       D2 = grid%ckwet(i2)


       select case(ip1)
         case(1)
           perm1 = perm_xx_loc_p(ng)
         case(2)
           perm1 = perm_yy_loc_p(ng)
         case(3)
           perm1 = perm_zz_loc_p(ng)
       end select
       
	   delxbc=0.D0
       select case(grid%ibndtyp(ibc))
          case(1)
		    delxbc=0.D0
          case(2)
          ! solve for pb from Darcy's law given qb /= 0
			  locpat%xxbc(:,nc)= xx((ng-1)*grid%ndof+1: ng*grid%ndof)
                delxbc=locpat%delx(1:grid%ndof,ng)
          case(3) 
          !    grid%xxbc(1,nc)=grid%pressurebc(2,ibc)
			  locpat%xxbc(2:grid%ndof,nc)= xx((ng-1)*grid%ndof+2: ng*grid%ndof)
			 ! grid%iphasebc(nc)=int(iphase_loc_p(ng))
			  delxbc(1)=0.D0
			  delxbc(2:grid%ndof)=locpat%delx(2:grid%ndof,ng)
		    end select

   
    iicap=  int(icap_loc_p(ng))     
       
	
	  call pri_var_trans_ims_ninc(locpat%xxbc(:,nc),temp,&
	     grid%scale,grid%nphase, &
      iicap,grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), & !use node's value
      grid%pcbetac(iicap),grid%pwrprm(iicap),&
      locpat%varbc(1:grid%size_var_use), ierr)
	
	
	
      call pri_var_trans_ims_winc(locpat%xxbc(:,nc), delxbc, temp,&
                   	   grid%scale,grid%nphase, &
					  iicap,grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
					  grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), & !use node's value
					  grid%pcbetac(iicap),grid%pwrprm(iicap),&
					  locpat%varbc(grid%size_var_use+1:(grid%ndof+1)*grid%size_var_use),ierr)
					  
	  !print *,' Mph Jaco BC terms: finish increment'
    do nvar=1,grid%ndof
   
    call IMSRes_FLBCCont(nc,locpat%areabc(nc), &
             locpat%varbc(nvar*grid%size_var_use+1:(nvar+1)*grid%size_var_use), &
			 var_loc_p((ng-1)*grid%size_var_node+nvar*grid%size_var_use+1:&
             (ng-1)*grid%size_var_node+nvar*grid%size_var_use+grid%size_var_use),&
			 porosity_loc_p(ng),tor_loc_p(ng),grid%sir(1:grid%nphase,iicap),locpat%distbc(nc),&
			 perm1,D2, locpat, grid, vv_darcy,Res)

    
	ResInc(m,1:grid%ndof,nvar) = ResInc(m,1:grid%ndof,nvar) - Res(1:grid%ndof)
   enddo

 enddo
   !print *,' Mph Jaco Finished BC terms'

 do n= 1, locpat%nlmax
   ra=0.D0
   ng = locpat%nL2G(n)
   na1= locpat%nG2N(ng) !PetSc Global index
   ! Remember, the matrix index starts from (0,0)
    p1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
   
   max_dev=0.D0
   do neq=1, grid%ndof
     do nvar=1, grid%ndof
       ra(neq,nvar)=ResInc(n,neq,nvar)/locpat%delx(nvar,ng)-ResOld_AR(n,neq)/locpat%delx(nvar,ng)
	   if(max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
	 enddo  	  
   enddo

   
   if(max_dev<1D-5)then
    !print *,'Mph Jaco max dev = ', max_dev
  endif
  
    
  ! Diagnal term
   do neq=1, grid%ndof
     do nvar=1, grid%ndof
        Jac_elem(n,0,neq,nvar)=ra(neq,nvar)
		Jac_Ind(n,0)=na1
     ! if(n==1)print *, neq, nvar, ra(neq,nvar)
        enddo
   enddo
        
 enddo

! -----------------------------contribution from transport----------------------


 ResInc=0.D0
 
 
  do nc = 1, locpat%nconn  ! For each interior connection...
    ra=0.D0
    m1 = locpat%nd1(nc) ! ghosted
    m2 = locpat%nd2(nc)

    n1 = locpat%nG2L(m1) ! = zero for ghost nodes
    n2 = locpat%nG2L(m2) ! Ghost to local mapping   
    na1= locpat%nG2N(m1)
	na2= locpat%nG2N(m2)
    !print *, grid%myrank,nc,m1,m2,n1,n2,na1,na2
    p1 =  (m1-1)*grid%ndof
    p2 =  (m2-1)*grid%ndof
   
    dd1 = locpat%dist1(nc)
    dd2 = locpat%dist2(nc)
    
    ip1 = locpat%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = locpat%iperm2(nc)
    
  
       i1 = ithrm_loc_p(m1)
       i2 = ithrm_loc_p(m2)
	    D1 = grid%ckwet(i1)
        D2 = grid%ckwet(i2)




    select case(ip1)
    case(1) 
       perm1 = perm_xx_loc_p(m1)
    case(2)
       perm1 = perm_yy_loc_p(m1)
    case(3)
       perm1 = perm_zz_loc_p(m1)
    end select
    
    select case(ip2)
    case(1) 
       perm2 = perm_xx_loc_p(m2)
    case(2)
       perm2 = perm_yy_loc_p(m2)
    case(3)
       perm2 = perm_zz_loc_p(m2)
    end select


    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
    iicap1=int(icap_loc_p(m1))
	iicap2=int(icap_loc_p(m2))
	
  ! do neq = 1, grid%ndof
    do nvar = 1, grid%ndof
    
		call IMSRes_FLCont(nc ,locpat%area(nc), &
         var_loc_p((m1-1)*grid%size_var_node+nvar*grid%size_var_use+1:&
         (m1-1)*grid%size_var_node+nvar*grid%size_var_use+grid%size_var_use),&
		 porosity_loc_p(m1),tor_loc_p(m1),grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1,&
		 var_loc_p((m2-1)*grid%size_var_node+1:(m2-1)*grid%size_var_node+grid%size_var_use),&
		 porosity_loc_p(m2),tor_loc_p(m2),grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,&
		 locpat, grid, vv_darcy,Res)

          ra(:,nvar)= Res(:)/locpat%delx(nvar,m1)-ResOld_FL(nc,:)/locpat%delx(nvar,m1)
	
			 
        call IMSRes_FLCont(nc ,locpat%area(nc), &
         var_loc_p((m1-1)*grid%size_var_node+1:(m1-1)*grid%size_var_node+grid%size_var_use),&
		 porosity_loc_p(m1),tor_loc_p(m1),grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1,&
		 var_loc_p((m2-1)*grid%size_var_node+nvar*grid%size_var_use+1:&
         (m2-1)*grid%size_var_node+nvar*grid%size_var_use+grid%size_var_use),&
		 porosity_loc_p(m2),tor_loc_p(m2),grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,&
		 locpat, grid, vv_darcy,Res)
 
          ra(:,nvar+grid%ndof)= Res(:)/locpat%delx(nvar,m2)-ResOld_FL(nc,:)/locpat%delx(nvar,m2)
	 
    enddo
  
   !   print *,' Mph Jaco Finished NC terms'
  
  ! enddo
   
 
      blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
  
	     do ii=0,grid%ndof-1
       do jj=0,grid%ndof-1
          if(n1>0) then
               blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
           endif
           if(n2>0) then
               blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
           endif
          enddo
   
		  do jj=grid%ndof,2*grid%ndof-1
              if(n1>0) then
                  blkmat12(ii+1,jj-grid%ndof+1) = blkmat12(ii+1,jj-grid%ndof+1) + ra(ii+1,jj+1)
              endif
              if(n2>0) then
                   blkmat22(ii+1,jj-grid%ndof+1) =  blkmat22(ii+1,jj-grid%ndof+1) - ra(ii+1,jj+1)
              endif
         enddo
    enddo
  
 ! then need determine index
   neighbor_ind1=0; neighbor_ind2=0;
   ndef=na1-na2

if(locpat%nlx>1)then
   if(ndef == -1)then
     neighbor_ind1=2; neighbor_ind2=1
   endif	 
   
   if(ndef == 1)then
     neighbor_ind1=1; neighbor_ind2=2
   endif	 
endif

if(locpat%nly>1)then
  if(ndef == -locpat%nlx)then
     neighbor_ind1=6; neighbor_ind2=5
	endif	 


   if(ndef == locpat%nlx)then
     neighbor_ind1=5; neighbor_ind2=6
	endif	 
endif


if(locpat%nlz>1)then
  if(ndef == locpat%nlxy)then
     neighbor_ind1=4; neighbor_ind2=3
	endif	 


   if(ndef == -locpat%nlxy)then
     neighbor_ind1=3; neighbor_ind2=4
	endif	 
endif
 
  
      if(n1>0)then
!	  call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
	  Jac_elem(n1,0,:,:)=Jac_elem(n1,0,:,:) + blkmat11(:,:)
      if(Jac_ind(n1,0)/=na1)print *, 'error @ index', n1,Jac_ind(n1,0),na1
	  endif
	  

	  if(n2>0)then
!	  call MatSetValuesBlocked(A,1,na2,1,na2,blkmat22,ADD_VALUES,ierr)
      Jac_elem(n2,0,:,:)=Jac_elem(n2,0,:,:) + blkmat22(:,:)
      if(Jac_ind(n2,0)/=na2)print *, 'error @ index', n2,Jac_ind(n2,0),na2
	  endif

	  
	  if(n1>0)then
	  !call MatSetValuesBlocked(A,1,na1,1,na2,blkmat12,ADD_VALUES,ierr)
	   Jac_elem(n1,neighbor_ind1,:,:)=Jac_elem(n1,neighbor_ind1,:,:) + blkmat12(:,:)
	   Jac_ind(n1,neighbor_ind1)=na2
	  endif
	  
	  if(n2>0)then
	   !  call MatSetValuesBlocked(A,1,na2,1,na1,blkmat21,ADD_VALUES,ierr)
       Jac_elem(n2,neighbor_ind2,:,:)=Jac_elem(n2,neighbor_ind2,:,:) + blkmat21(:,:)
	   Jac_ind(n2,neighbor_ind2)=na1
	  endif 
 
 !  print *, nc,n1,n2,na1,na2, neighbor_ind1,neighbor_ind2,Jac_ind(na1,:),Jac_ind(na2,:)  
!print *,'accum r',ra(1:5,1:8)   
 !print *,'devq:',nc,q,dphi,devq(3,:)
  end do
  !print *,' IMS Jaco Finished Two node terms'
 
  ! nullify(Resold_AR, Resold_FL)
  
  
 end subroutine IMSJacobin_patch

end module IMS_patch_module


