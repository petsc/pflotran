
!#include "pflowgrid_ptran_init.F90"
module pflotran_couple_module
use pflow_gridtype_module
  
implicit none
private


#include "definitions.h"

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"


public pflowGrid_ptran_init, pflotranGrid_interp, ptran_bc_reassign


contains
  subroutine pflowGrid_ptran_init(grid, ierr)

  use ptran_global_module
  use trdynmem_module
  
  implicit none

  type(pflowGrid), intent(inout) :: grid
  
  integer, intent(in) :: ierr

  real*8, pointer :: por_pflow_p(:), por0_pflow_p(:),por_ptran_p(:), & !volume_p(:), 
                     temp_pflow_p(:), temp_ptran_p(:), &
                     sat_pflow_p(:), sat_ptran_p(:)
  real*8 :: afac
  integer*4 :: jn,n,nr
  
  grid%using_pflowGrid = PETSC_TRUE
  using_pflowGrid = PETSC_TRUE

! assignments by value
!  nx = grid%nx; ny = grid%ny; nz = grid%nz
!  npx = grid%npx; npy = grid%npy; npz = grid%npz
  
  nconn = grid%nconn
  nconnbc=grid%nconnbc

  allocate(vl(nconn))
  allocate(vlbc(nconnbc))
  if (iphase == 2) then
    allocate(vg(nconn))
    allocate(vgbc(nconnbc))
  endif
  
! pointer assignments
  nL2G  => grid%nL2G
  nG2L  => grid%nG2L
  NL2A  => grid%nL2A
  nd1   => grid%nd1
  nd2   => grid%nd2
  dist1 => grid%dist1
  dist2 => grid%dist2
  area  => grid%area
  nblkbc = grid%nblkbc
  distbc => grid%distbc
  areabc => grid%areabc

  ! nblkbc need to be changed
  mblkbc => grid%mblkbc
  ibconn => grid%ibconn

! set plotting times
  kplot = grid%kplot
  do n = 1, grid%kplot
    tplot(n) = grid%tplot(n)
  enddo
  
! allocate(vb(grid%nlmax))
! call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, por_pflow_p, ierr)
  call VecGetArrayF90(grid%porosity0, por0_pflow_p, ierr)
  call VecGetArrayF90(porosity, por_ptran_p, ierr)
  call VecGetArrayF90(grid%temp, temp_pflow_p, ierr)
  call VecGetArrayF90(temp, temp_ptran_p, ierr)
  call VecGetArrayF90(grid%sat, sat_pflow_p, ierr)
  call VecGetArrayF90(sat, sat_ptran_p, ierr)
  !call VecGetArrayF90(phik, phik_ptran_p, ierr)
  do n = 1, grid%nlmax
    jn = 2+(n-1)*grid%nphase
!   vb(n) = volume_p(n)
    temp_ptran_p(n) = temp_pflow_p(n)
    por_ptran_p(n) = por_pflow_p(n)
    afac = (1.d0-por_pflow_p(n))/(1.d0-por0_pflow_p(n))
    do nr = 1, nkin
      phik_p(nr+(n-1)*nkin) = afac*phik_p(nr+(n-1)*nkin)
!     if (phik_p(nr+(n-1)*nkin) > 1.d-8) &
!     surf_p(nr+(n-1)*nkin) = surf_p(nr+(n-1)*nkin)*phik_p(nr+(n-1)*nkin)
      
      surf_p(nr+(n-1)*nkin) = 1.d0
      
!     print *,'couple: ',n,nr,phik_p(nr+(n-1)*nkin),surf_p(nr+(n-1)*nkin)
    enddo
!   phik_p((n-1)*nkin+1: n*nkin)=phik_p((n-1)*nkin+1: n*nkin)&
!          *(1.D0-por_pflow_p(n)**2.D0/por0_pflow_p(n))/(1.D0-por_pflow_p(n))
!   print *, 'couple:', por0_pflow_p(n),por_pflow_p(n), phik_p((n-1)*nkin+1: n*nkin)

    sat_ptran_p(n) = 1.d0-sat_pflow_p(jn)

!   print *,'pflowGrid_ptran_init: ',n,sat_ptran_p(n),sat_pflow_p(n),por_pflow_p(n)
  enddo
! call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, por_pflow_p, ierr)
  call VecRestoreArrayF90(grid%porosity0, por0_pflow_p, ierr)
  !call VecRestoreArrayF90(phik, phik_ptran_p, ierr)
  call VecRestoreArrayF90(porosity, por_ptran_p, ierr)
  call VecRestoreArrayF90(grid%temp, temp_pflow_p, ierr)
  call VecRestoreArrayF90(temp, temp_ptran_p, ierr)
  call VecRestoreArrayF90(grid%sat, sat_pflow_p, ierr)
  call VecRestoreArrayF90(sat, sat_ptran_p, ierr)
  
  surf0_p(:) = surf_p(:)
  phik0_P(:) = phik_p(:)
  !call VecCopy(phik, phik0, ierr)
  call VecCopy(grid%porosity,grid%porosity0,ierr)

  end subroutine pflowGrid_ptran_init

!======================================================================

  subroutine pflotranGrid_interp (grid, tflow1, tflow2)

  use ptran_global_module
  use trdynmem_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
    
  integer*4 :: nc, m, jm
  integer :: ierr, iiphas 
  real*8 :: f1, f2, tflow1, tflow2

  real*8, pointer :: xx_p(:), yy_p(:),iphase_p(:),iphase_old_p(:),&
                     temp_p(:), ttemp_p(:), ptran_temp_p(:), &
                     press_p(:), ppress_p(:), ptran_press_p(:),&
                     sat_p(:),ssat_p(:), ptran_sat_p(:),den_co2_p(:), xphi_co2_p(:)

! interpolate field variables at time t
  f1 = (t - tflow1)/(tflow2 - tflow1)
! f2 = (tflow2 - t)/(tflow2 - tflow1)
  f2 = 1.d0 - f1

  call VecGetArrayF90(grid%temp,temp_p,ierr)
  call VecGetArrayF90(grid%ttemp,ttemp_p,ierr)
  call VecGetArrayF90(temp,ptran_temp_p,ierr)

  call VecGetArrayF90(grid%pressure,press_p,ierr)
  call VecGetArrayF90(grid%ppressure,ppress_p,ierr)
  call VecGetArrayF90(press,ptran_press_p,ierr)

  call VecGetArrayF90(grid%sat,sat_p,ierr)
  call VecGetArrayF90(grid%ssat,ssat_p,ierr)
  call VecGetArrayF90(ssat,ptran_sat_p,ierr)

  if (grid%use_2ph == PETSC_TRUE) then
!   call VecGetArrayF90(grid%xx_loc, xx_p, ierr)
    call VecGetArrayF90(grid%xx, xx_p, ierr)
    do m = 1, grid%nlmax
!     n = nL2G(m)
!     ppress_p(m) = xx_p(1+(n-1)*grid%ndof) ! ppress and press defined for nphase
!     ttemp_p(m) = xx_p(2+(n-1)*grid%ndof)
      jm = 2+(m-1)*grid%nphase 
      ppress_p(m) = xx_p(1+(m-1)*grid%ndof) ! ppress and press defined for nphase
      ttemp_p(m) = xx_p(2+(m-1)*grid%ndof)
      ssat_p(jm) = xx_p(4+(m-1)*grid%ndof)
	  enddo
!   call VecRestoreArrayF90(grid%xx_loc, xx_p, ierr)
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  endif

  if (grid%use_mph == PETSC_TRUE) then
!   call VecGetArrayF90(grid%xx_loc, xx_p, ierr)
    call VecGetArrayF90(grid%xx, xx_p, ierr)
    call VecGetArrayF90(grid%yy, yy_p, ierr)
    call VecGetArrayF90(grid%iphas, iphase_p, ierr)
    call VecGetArrayF90(grid%iphas_old, iphase_old_p, ierr)	
    do m = 1, grid%nlmax
!     n = nL2G(m)
!     ppress_p(m) = xx_p(1+(n-1)*grid%ndof) ! ppress and press defined for nphase
!     ttemp_p(m) = xx_p(2+(n-1)*grid%ndof)
      jm = 2+(m-1)*grid%nphase 
      ppress_p(m) = xx_p(1+(m-1)*grid%ndof) ! ppress and press defined for nphase
      press_p(m) = yy_p(1+(m-1)*grid%ndof)
      ttemp_p(m) = xx_p(2+(m-1)*grid%ndof)
      temp_p(m) = yy_p(2+(m-1)*grid%ndof)
      
      iiphas=iphase_p(m)
      select case(iiphas)
        case(1)
          ssat_p(jm) = 0.D0
        case(2)
          ssat_p(jm) = 1.D0
        case(3)
          ssat_p(jm) = xx_p(3+(m-1)*grid%ndof)    	  
      end select
       
      iiphas=iphase_old_p(m)
      select case(iiphas)
        case(1)
          sat_p(jm) = 0.D0
        case(2)
          sat_p(jm) = 1.D0 
        case(3)
          sat_p(jm) = yy_p(3+(m-1)*grid%ndof)    	  
      end select
    enddo
!   call VecRestoreArrayF90(grid%xx_loc, xx_p, ierr)
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
    call VecRestoreArrayF90(grid%yy, yy_p, ierr)
    call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
    call VecRestoreArrayF90(grid%iphas_old, iphase_old_p, ierr)	
  endif

  call VecGetArrayF90(den_co2, den_co2_p, ierr)
  call VecGetArrayF90(xphi_co2, xphi_co2_p, ierr)
  do m = 1, grid%nlmax
    jm = 2+(m-1)*grid%nphase
    ptran_temp_p(m) = f1 * ttemp_p(m) + f2 * temp_p(m)
    ptran_press_p(m) = f1 * ppress_p(m) + f2 * press_p(m)
    
    ptran_sat_p(m) = f1 *(1.D0-ssat_p(jm)) + f2 * (1.D0-sat_p(jm))
    
    if (ptran_sat_p(m) < slcutoff) ptran_sat_p(m)= 0.9D0 * slcutoff
    
    den_co2_p(m) =  f1 * grid%dden_co2(m) + f2 * grid%den_co2(m)
    xphi_co2_p(m) = f1 * grid%xxphi_co2(m) + f2 * grid%xphi_co2(m)
    
!  if (m == 25) &
!  print *,'pflotranGrid_interp: ',grid%myrank,m,ptran_temp_p(m), &
!  ttemp_p(m),temp_p(m),ptran_press_p(m),ppress_p(m),press_p(m), &
!  ptran_sat_p(m),ssat_p(jm),sat_p(jm),f1,f2,xphi_co2(m),grid%xxphi_co2(m),grid%xphi_co2(m)

!   if (m < 5) &
!   print *,'pflotranGrid_interp: ',grid%myrank,m,ptran_temp_p(m),ptran_press_p(m),&
!   ptran_sat_p(m),f1,f2,xphi_co2(m),grid%xxphi_co2(m),grid%xphi_co2(m)
  enddo

  call VecRestoreArrayF90(grid%temp,temp_p,ierr)
  call VecRestoreArrayF90(grid%ttemp,ttemp_p,ierr)
  call VecRestoreArrayF90(temp,ptran_temp_p,ierr)
  
  call VecRestoreArrayF90(grid%pressure,press_p,ierr)
  call VecRestoreArrayF90(grid%ppressure,ppress_p,ierr)
  call VecRestoreArrayF90(press,ptran_press_p,ierr)
  
  call VecRestoreArrayF90(grid%sat,sat_p,ierr)
  call VecRestoreArrayF90(grid%ssat,ssat_p,ierr)
  call VecRestoreArrayF90(ssat,ptran_sat_p,ierr)
  call VecRestoreArrayF90(den_co2, den_co2_p, ierr)
  call VecRestoreArrayF90(xphi_co2, xphi_co2_p, ierr)


! print *,'pflowgrid_ptran_init: ',f1,f2
#if 0
  if (f2 > 0.d0) then
    vl   = f1 * grid%vvl_loc + f2 * grid%vl_loc
    vlbc = f1 * grid%vvlbc   + f2 * grid%vlbc
    if(iphase==2)then
      vg   =  f1 * grid%vvg_loc + f2 * grid%vg_loc
      vgbc =  f1 * grid%vvgbc   + f2 * grid%vgbc
    endif
  else ! first step
    vl   = grid%vvl_loc
    vlbc = grid%vvlbc
    if(iphase==2)then
      vg   = grid%vvg_loc
      vgbc = grid%vvgbc
    endif
  endif
#endif

  call DAGlobalToLocalBegin(grid%da_1_dof, xphi_co2, INSERT_VALUES, &
                            xphi_co2_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, xphi_co2, INSERT_VALUES, &
                            xphi_co2_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, den_co2, INSERT_VALUES, &
                            den_co2_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, den_co2, INSERT_VALUES, &
                            den_co2_loc, ierr)

  if (f2 > 0.d0) then
    do nc = 1, grid%nconn
      vl(nc) = f1 * grid%vvl_loc(nc) + f2 * grid%vl_loc(nc)
      if(iphase==2) vg(nc) = f1 * grid%vvg_loc(nc) + f2 * grid%vg_loc(nc)
      
!     print *,'interp-l: ',nc,vl(nc),f1,f2,grid%vvl_loc(nc),grid%vl_loc(nc)
!     print *,'interp-g: ',nc,vg(nc),f1,f2,grid%vvg_loc(nc),grid%vg_loc(nc)
    enddo
    do nc = 1, grid%nconnbc
      vlbc(nc) = f1 * grid%vvlbc(nc) + f2 * grid%vlbc(nc)
	  if(iphase==2)vgbc(nc) = f1 * grid%vvgbc(nc) + f2 * grid%vgbc(nc)
    enddo
  else ! first step or f2 = 0, f1 = 1
    do nc = 1, grid%nconn
      vl(nc) = grid%vvl_loc(nc)
      if(iphase==2) vg(nc) = grid%vvg_loc(nc)
      
!     print *,'interp-l: ',nc,vl(nc),f1,f2,grid%vvl_loc(nc),grid%vl_loc(nc)
!     print *,'interp-g: ',nc,vg(nc),f1,f2,grid%vvg_loc(nc),grid%vg_loc(nc)
    enddo
    do nc = 1, grid%nconnbc
      vlbc(nc) = grid%vvlbc(nc)
      if(iphase==2) vgbc(nc) = grid%vvgbc(nc)
    enddo
  endif
! print *, f1,f2,vlbc,grid%vvlbc,grid%vlbc
  end subroutine pflotranGrid_interp
  

!========================================================================================
  
  subroutine ptran_bc_reassign(grid)

      use ptran_global_module
      use trdynmem_module
      use water_eos_module
      use ptran_speciation_module

      implicit none
	  type(pflowGrid), intent(in) :: grid

	  integer iibndtyp(nrgbcmx)

      real*8 :: cloc(ncmx),cxloc(ncxmx),cloctot(ncmx),pgasloc(ngmx), &
      cecloc(nexmx),xexloc(nexmx), &
      ccsorplc(nscxmx),csorplc(nscxmx),csorpflc(nscxmx), &
      sitdnloc(nscxmx), &
      dgamdi(ncmx),gamloc(ncmx),gamxloc(ncxmx)

      real*8 :: fach2o,rho1,sum,u1
      real*8 ctot_ori(ncmx,nrgmx), prefbc(nrgmx)
	  integer  ibndtyp_ori(nrgbcmx), iface_ori(nrgbcmx), itype_ori(ncmx,nrgmx)
	  character(len=namlen) :: ncon_ori(ncmx,nrgmx)
	   
      integer :: i,idif=0,iireg,iisrc,iter,j,jj,k,l,m,nss
      integer iibc, ptran_ibcreg1,ptran_ibcreg2

      character :: psinam*7

	   
	   print *,'pflow bc number:',grid%nblkbc
!	   do ibc=1,grid%nblkbc
!	    print*,'pflotran_couple: ', ibc, grid%xxbc0(1:2,ibc)
!	   enddo	
	   

      ibc = 0
        
	  
	  
	  
	  
    if(ibcreg(1).eq.0) goto 110 ! source/sink
      
	  
	   iibc=0;ibc=0
       ptran_ibcreg1=ibcreg(1)
	   ptran_ibcreg2=ibcreg(2)
	   ibcreg(2)=ibcreg(1)+grid%nblkbc -1
	   ibndtyp_ori(:)=ibndtyp(:) 
	   ctot_ori(:,:)=ctot(:,:)
	   itype_ori(:,:)=itype(:,:)
	   iface_ori(:)=iface(:)
	   ncon_ori(:,:)=ncon(:,:)
	   
	    do ireg = ptran_ibcreg1, ptran_ibcreg2
	       ibc=ibc+1 !now ibc+ptran_ibcreg1=ireg
		   
!          print *, ireg,ibc, iface_ori(ibc)
!          print *,'coup_ori:ibc=', ibc, ireg,iibc, iface_ori(ibc), ibndtyp_ori(ibc),&
!		             ctot_ori(9:10,ireg)
		   
		   if(iface_ori(ibc)==1 .or. iface_ori(ibc)==2)then
		    
			
			
			
			do i=1, grid%nz
		       iibc=iibc+1
               iface(iibc)=iface_ori(ibc)
			   ibndtyp(iibc)=ibndtyp_ori(ibc)
!              print *,'pflotran_couple: ', iibc,iface_ori(ibc),iface(iibc) 
!              print *, iibc,grid%xxbc0(1:2,iibc) 
		       tempbc(iibc)=grid%xxbc0(2,iibc) + tkelvin
			   prefbc(iibc)=grid%xxbc0(1,iibc)
			   !print *,iibc, iface_ori(ireg), ibndtyp_ori(ibc)
			   itype(:,ptran_ibcreg1+iibc-1)=itype_ori(:,ireg)
			   ctot(:,ptran_ibcreg1+iibc-1)=ctot_ori(:,ireg)
			   ncon(:,ptran_ibcreg1+iibc-1)=ncon_ori(:,ireg)
             enddo
            endif
		  
				   
		   if(iface_ori(ibc)==3 .or. iface_ori(ibc)==4)then
		       iibc=iibc+1
               iface(iibc)=iface_ori(ibc)
			   ibndtyp(iibc)=ibndtyp_ori(ibc)
		         tempbc(iibc)=grid%xxbc0(2,iibc) + tkelvin
			   prefbc(iibc)=grid%xxbc0(1,iibc)
			   itype(:,ptran_ibcreg1+iibc-1)=itype_ori(:,ireg)
			   ctot(:,ptran_ibcreg1+iibc-1)=ctot_ori(:,ireg)
			   ncon(:,ptran_ibcreg1+iibc-1)=ncon_ori(:,ireg)
!             print *,'pflotran_couple: ', iibc,ptran_ibcreg1+iibc, ctot(:, ptran_ibcreg1+iibc-1)
			 endif
      
           if(iface_ori(ibc)==5 .or. iface_ori(ibc)==6)then
		    
			! front & back not working even in pflow part!!!
			! only type 2 allowed i.e. ibndtyp==2 ! 
		       iibc=iibc+1
               iface(iibc)=iface_ori(ibc)
			   ibndtyp(iibc)=ibndtyp_ori(ibc)
		       tempbc(iibc)=grid%xxbc0(2,iibc) + tkelvin
			   prefbc(iibc)=grid%xxbc0(1,iibc)
			   itype(:,ptran_ibcreg1+iibc-1)=itype_ori(:,ireg)
			   ctot(:,ptran_ibcreg1+iibc-1)=ctot_ori(:,ireg)
			   ncon(:,ptran_ibcreg1+iibc-1)=ncon_ori(:,ireg)
          
           endif
      enddo
    
     				  		      



     


      do ibc = 1, grid%nblkbc
	            ireg=ibc + ptran_ibcreg1 -1
				call density(tempbc(iibc), prefbc(iibc),rho1)
			   	rho(1)=1.D0!rho1
				 	        
			    call trspeciate (cloc,cxloc,cloctot,pgasloc, &
                 cecloc,xexloc,csorplc,ccsorplc,csorpflc,sitdnloc,alogpf, &
                 gamloc,gamxloc,dgamdi,rho,tempbc(ibc),iter)

                               fach2o = one
                   if (molal.eq.0) fach2o = one/(wh2o*cloc(jh2o))


                    psinam = ' psibnd'
                     do j = 1, ncomp
            !         psibnd(j,ibc) = cloctot(j) * rho1
                      psibnd(j,ibc) = cloctot(j)
                      ccbnd (j,ibc) = cloc(j)
                     if (myrank==0) &
                     write(iunit2,1020) nam(j),psinam,j,ibc,psibnd(j,ibc)
                     enddo

                   psinam = ' xexbnd'
                  do m = 1, nexsolid
                      do k = nsitex1(m), nsitex2(m)
                           do jj = nex1(k), nex2(k)
                               j = jex(jj)
                               xexbnd(jj,ibc) = xexloc(jj)
                               if (myrank==0) &
                                write(iunit2,1020) nam(j),psinam,jj,ibc,xexbnd(jj,ibc)
                           enddo
                        enddo
                   enddo

				if (nsrfmin .gt. 0) then
				  psinam = ' psorpbnd'
				  do j = 1, ncomp
					sum = zero
					do m = nsrfmin-ncolsrf+1,nsrfmin ! sum over colloids only !
					  do k = nsite1(m), nsite2(m)
						do i = nsorp1(k), nsorp2(k)
						  sum = sum + ssorp(j,i)*ccsorplc(i)
						enddo
					  enddo
					enddo
					psorpbnd(j,ibc) = sum
					if (myrank==0) &
					write(iunit2,1020) nam(j),psinam,j,ibc,sum
				  enddo
				endif

        if (iphase.eq.2 .or. iphase.eq.0) then !  conc for the gas phase
          u1 = one/(rgasjj*tk)
          psinam = 'psigbnd'
          do j = 1, ncomp
            sum = zero
            do l = 1, ngas
              sum = sum + sgas(j,l)*pgasloc(l)
            enddo
            psigbnd(j,ibc) = sum*u1                           
            if (myrank==0) &
            write(iunit2,1020) nam(j),psinam,j,ibc,psigbnd(j,ibc)
          enddo
          do i = 1, ngas
            pgasbnd(i,ibc) = pgasloc(i)*u1                        
          enddo
        endif
!       print *,'pflotran_couple: ',ibc, psibnd(9:10,ibc)
  
    enddo
  100 continue
  
  
    
   deallocate(grid%xxbc0)
   deallocate(grid%iphasebc0)
  
  110 return
1010 format (5x,'--> boundary condition:  type =',i2, &
      '  temp = ',1pg10.3,'  iter =',i4)   
1020 format(10x,'component = ',a12,3x,a7,'(',i2,i3,') =',1p10e12.4)


  end subroutine ptran_bc_reassign

  end module pflotran_couple_module
	