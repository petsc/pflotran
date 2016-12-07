! introduced grid variables: e_total :: 1 dof
! translator_module                           c_total :: grid%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (grid%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

#include "include/finclude/petscsnes.h"
  use petscsnes
  module IMS_module
  use pflow_gridtype_module
 ! use pflow_var_module

private 

  public IMSResidual, IMSJacobin, pflow_IMS_initaccum, &
         pflow_update_IMS,pflow_IMS_initadj, pflow_IMS_timecut,&
         pflow_IMS_setupini, IMS_Update_Reason ,&
		 pflow_ims_step_maxchange



  contains


 subroutine pflow_ims_massbal(grid, locpat)
 use pflow_gridtype_module
  implicit none
  type(pflowGrid) :: grid 
  type(pflow_localpatch_info), intent(inout) :: locpat
   
  integer :: ierr,icall
  integer :: n,n0,nc,np,n2p,n2p0,ng
  real*8 x,y,z,nzm,nzm0, nxc,nxc0,c0, c00,nyc,nyc0,nzc,nzc0,nsm,nsm0,sm 
  integer :: index
     
  PetscScalar, pointer ::  porosity_p(:), volume_p(:),  var_p(:)
                           
   
  real*8 ::  pvol,sum
  real*8, pointer ::  den(:),sat(:)
 
  real*8 :: tot(0:grid%nphase), tot0(0:grid%nphase)  
  data icall/0/

  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  var_p => locpat%var
 
  tot=0.D0
  n2p=0
  nxc=0.; nyc=0.; nzc=0.D0; nzm=grid%z(grid%nmax); sm=0.; nsm=nzm; c0=0.D0

  do n = 1,locpat%nlmax
    n0=(n-1)* grid%ndof
	ng=locpat%nL2G(n)
    index=(ng-1)*grid%size_var_node
    den=>var_p(index+3+grid%nphase: index+2+2*grid%nphase)
    sat=>var_p(index+2+1:index+2+grid%nphase)
    
   
    pvol=volume_p(n)*porosity_p(n)		     
	 n2p=n2p+1
     x=grid%x(locpat%nL2A(n)+1)
	 y=grid%y(locpat%nL2A(n)+1)
	 z=grid%z(locpat%nL2A(n)+1)
	 
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
  
  nullify(var_p)
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
   &" gas:",1p, e13.6, " tot:", 1p, e13.6, " [kg]")') &
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
    write(13,'(1p5e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,&
 	tot(1),tot(2),tot(1)+tot(2)!,real(n2p),nzc,nzm,nsm 
  endif    
  
  
  
  
 end subroutine pflow_ims_massbal


 subroutine pflow_ims_step_maxchange(grid)
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

 end  subroutine pflow_ims_step_maxchange
 
  
 
   

 subroutine pflow_IMS_timecut(grid, locpat)
 
  implicit none
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
  
 
  PetscScalar, pointer :: xx_p(:),yy_p(:)!,var_p(:),iphase_p(:)
  integer :: n,n0,re,ierr
  !integer re0, ierr, index, iipha
  !real*8, pointer :: sat(:),xmol(:)

   call VecGetArrayF90(grid%xx, xx_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr)
 ! call VecGetArrayF90(grid%var, var_p, ierr); 
 ! call VecGetArrayF90(grid%iphas, iphase_p, ierr); 

  do n=1, locpat%nlmax
   n0=(n-1)*grid%ndof
   do re = 1, grid%ndof
   !xx_p(n0+re)= 0.5D0 * xx_p(n0+re) +.5D0 *yy_p(n0+re)
   xx_p(n0+re)= yy_p(n0+re)
   enddo
   enddo 
   call VecRestoreArrayF90(grid%xx, xx_p, ierr) 
    call VecRestoreArrayF90(grid%yy, yy_p, ierr)
	
	!call VecCopy(grid%xx,grid%yy,ierr)
	!call pflow_IMS_initaccum(grid)
 
  end subroutine pflow_IMS_timecut
  

 subroutine pflow_IMS_setupini(grid, locpat)
  implicit none
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
  
interface
   subroutine pims_vecgetarrayf90(grid, patch, vec, f90ptr, ierr)
     use pflow_gridtype_module
     implicit none
#include "include/finclude/petsc.h"

     type(pflowGrid), intent(inout) :: grid
     type(pflow_localpatch_info) :: patch
     Vec :: vec
     PetscScalar, dimension(:), pointer :: f90ptr
     integer :: ierr
   end subroutine pims_vecgetarrayf90
end interface

  PetscScalar, pointer :: xx_p(:)
  integer iln,na,nx,ny,nz,ir,ierr
  
   grid%size_var_use = 2 + 4*grid%nphase
    grid%size_var_node = (grid%ndof + 1) * grid%size_var_use


   allocate(locpat%var(locpat%ngmax * grid%size_var_node))	
   allocate(locpat%delx(grid%ndof,locpat%ngmax))
   locpat%delx=0.D0
   
  call pims_vecgetarrayf90(grid, locpat, grid%xx, xx_p, ierr); CHKERRQ(ierr)
  
  do iln=1, locpat%nlmax
    if(grid%Samrai_drive==PETSC_TRUE) then
         nz = int(iln/locpat%nlxy) + 1
         ny = int(mod(iln,locpat%nlxy)/locpat%nlx) + 1
         nx = mod(mod(iln,locpat%nlxy),locpat%nlx) + 1  
     else
     na = locpat%nL2A(iln)
    !compute i,j,k indices from na: note-na starts at 0
      nz = na/grid%nxy + 1
      ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
      nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
    endif
!   print *,'pflow_IMS_resjac: ',na,nx,ny,nz
    
    do ir = 1,grid%iregini
      if ((nz>=grid%k1ini(ir)) .and. (nz<=grid%k2ini(ir)) .and.&
          (ny>=grid%j1ini(ir)) .and. (ny<=grid%j2ini(ir)) .and.&
          (nx>= grid%i1ini(ir)) .and. (nx<=grid%i2ini(ir)))then
                 xx_p(1+(iln-1)*grid%ndof: iln*grid%ndof)=grid%xx_ini(:,ir)
		!	print *,'Setup ini ::', na,xx_p(1+(iln-1)*grid%ndof: iln*grid%ndof) 		   
                !exit
           endif
	enddo 
 enddo
 
 if(grid%Samrai_drive==PETSC_FALSE) then						
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
 endif

end  subroutine pflow_IMS_setupini
  


! Subroutine for solution checking 
 subroutine IMS_Update_Reason(reason,grid, locpat)
  
  implicit none
 
  integer, intent(out):: reason
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
  PetscScalar, pointer :: xx_p(:),var_p(:), yy_p(:),r_p(:)
  integer :: n,n0,re, i
  integer re0, ierr, index, iipha
  real*8 rmax(grid%ndof)

  !call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  re=1
 ! call SNESComputeFunction(grid%snes,grid%xx,grid%r,ierr)
 ! do n=1,grid%ndof
 !  call VecStrideNorm(grid%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if(re>0)then
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr)
  
  do n = 1,locpat%nlmax
     n0=(n-1)* grid%ndof
  	 do i=2, grid%ndof
   	   if(xx_p(n0 + i) > 1.0D0)then
	    re=0; exit
!		goto 111
       endif
     enddo
	if(re<=0) exit
  enddo 	  		

 ! print *, 'update reason: ',grid%myrank,grid%nlmax,n,re

  !   call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
	 !print *,' update reason ba MPI', ierr
  if(re<=0) print *,'Sat out of Region at: ',n,xx_p(n0+1:n0+grid%ndof)
    call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(grid%yy, yy_p, ierr)
 endif

  
  if(grid%commsize >1)then
    call MPI_ALLREDUCE(re, re0,1, MPI_INTEGER,MPI_SUM, &
    PETSC_COMM_WORLD,ierr)
	!print *,' update reason re'
    !call MPI_BCAST(re0,1, MPI_INTEGER, 0,PETSC_COMM_WORLD,ierr)
	!print *,' update reason ca'
    if(re0<grid%commsize) re=0
  endif
  reason=re

  end subroutine IMS_Update_Reason


! Residual called by SNESSolve
  subroutine IMSResidual(snes,xx,r,grid, ierr)
    use translator_IMS_module
    use IMS_patch_module
    implicit none
 
    SNES, intent(in) :: snes
    Vec, intent(inout) :: xx
    Vec, intent(out) :: r
    type(pflowGrid), intent(inout) :: grid
    type(pflow_localpatch_info), pointer:: locpat
    PetscScalar, pointer :: xx_p(:), accum_p(:),r_p(:)
     
   integer :: ierr, n,i
 
  !grid%vvlbc=0.D0
  !grid%vvgbc=0.D0
  !grid%vvl_loc=0.D0
  !grid%vvg_loc=0.D0

  ! print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

  !if(grid%dt>5D5)
   locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr
   !print *,'Res ims :: in res', locpat%nlmax

	call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
     do n = 1, locpat%nlmax
	    do i=2,grid%nphase
		   if(xx_p((n-1)*grid%ndof+i) <0.D0)xx_p((n-1)*grid%ndof+i) =1.D-16
		   if(xx_p((n-1)*grid%ndof+i) >1.D0)xx_p((n-1)*grid%ndof+i) =1.D0 -1D-16
        enddo
      
     enddo
    call VecRestoreArrayF90(xx, xx_p, ierr)
	
  !print *,' IMS_RES : Begin scatter'
	
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
   
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_xx, &
                            INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_xx, &
                          INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_yy, &
                            INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_yy, &
                          INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_zz, &
                            INSERT_VALUES, grid%perm_zz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_zz, &
                          INSERT_VALUES, grid%perm_zz_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%ithrm, &
                            INSERT_VALUES, grid%ithrm_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%ithrm, &
                          INSERT_VALUES, grid%ithrm_loc, ierr)
   call DAGlobalToLocalBegin(grid%da_1_dof, grid%icap, &
                            INSERT_VALUES, grid%icap_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%icap, &
                          INSERT_VALUES, grid%icap_loc, ierr)


 !print *,' IMS_RES : End scatter'
! End distribute data 
! now assign access pointer to local variables
  call VecGetArrayF90(grid%xx_loc, locpat%xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
! call VecGetArrayF90(grid%yy, yy_p, ierr)
 

  ! notice:: here we assume porosity is constant
 
  call VecGetArrayF90(grid%yy,locpat%yy_p,ierr)
  call VecGetArrayF90(grid%porosity_loc, locpat%porosity_loc_p, ierr)
  call VecGetArrayF90(grid%tor_loc, locpat%tor_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, locpat%perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, locpat%perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, locpat%perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, locpat%volume_p, ierr)
  call VecGetArrayF90(grid%ithrm_loc, locpat%ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, locpat%icap_loc_p, ierr)
  call VecGetArrayF90(grid%vl, locpat%vl_p, ierr)
  
 ! print *,' IMS_res : Finished scattering non deriv'


! Evaluate residual on a patch
!  In petsc version, only one patch on a processor
!   all processor will call IMSResidual_patch, not collective    
  call IMSResidual_patch(locpat%xx_loc_p,r_p,accum_p,locpat,grid,ierr)


  
 call VecRestoreArrayF90(r, r_p, ierr)
 call VecRestoreArrayF90(grid%accum, accum_p, ierr)
 call VecRestoreArrayF90(grid%yy, locpat%yy_p, ierr)
 call VecRestoreArrayF90(grid%xx_loc, locpat%xx_loc_p, ierr)
 call VecRestoreArrayF90(grid%porosity_loc, locpat%porosity_loc_p, ierr)
 call VecRestoreArrayF90(grid%tor_loc, locpat%tor_loc_p, ierr)
 call VecRestoreArrayF90(grid%perm_xx_loc, locpat%perm_xx_loc_p, ierr)
 call VecRestoreArrayF90(grid%perm_yy_loc, locpat%perm_yy_loc_p, ierr)
 call VecRestoreArrayF90(grid%perm_zz_loc, locpat%perm_zz_loc_p, ierr)
 call VecRestoreArrayF90(grid%volume, locpat%volume_p, ierr)
 call VecRestoreArrayF90(grid%ithrm_loc, locpat%ithrm_loc_p, ierr)
 call VecRestoreArrayF90(grid%icap_loc, locpat%icap_loc_p, ierr)
 call VecRestoreArrayF90(grid%vl, locpat%vl_p, ierr)
 

!print *,'Residual ::...........'; call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)
! print *,'finished IMSResidual'
 end subroutine IMSResidual
                
! --------------------------------------------------------------------- 

  subroutine IMSJacobin(snes,xx,A,B,flag,grid, ierr)
       
    use translator_IMS_module
    use IMS_patch_module
    
    implicit none

    SNES, intent(in) :: snes
    Vec, intent(in) :: xx
    Mat, intent(out) :: A, B
    type(pflowGrid), intent(inout) :: grid
    type(pflow_localpatch_info), pointer :: locpat
   ! integer, intent(inout) :: flag
    MatStructure flag

    integer :: ierr, n,ng, na1,na2, ind,neq,nvar
   
    real*8 :: blkmat11(1:grid%ndof,1:grid%ndof)

   ! matrix block, 0 to 6 should be 3D non-zero entrys.!
	!               0 stands for diagnal element
	!               1 left.  2. right  3 up.  4 down.  5 front. 6 back.
	!               This is determined by the connection setup, does not really matter   

  	real*8, pointer :: Jac_elem(:,:,:,:)

	! non-zero matrix cloumn index, diagnal will be the same as row index
	! for petsc version, it should be global index 
	integer,pointer :: Jac_Ind(:,:)

    
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
 flag = SAME_NONZERO_PATTERN
 locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr


	allocate(Jac_elem(1:locpat%nlmax,0:6,1:grid%ndof,1:grid%ndof))
	allocate(Jac_Ind(1:locpat%nlmax,0:6))
 
 
  !print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
 ! call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
 !                           grid%xx_loc, ierr)
 ! call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
 !                         grid%xx_loc, ierr)
 
 !call DAGlobalToLocalBegin(grid%da_1_dof, grid%iphas, &
 !                          INSERT_VALUES, grid%iphas_loc, ierr)
 !call DAGlobalToLocalEnd(grid%da_1_dof, grid%iphas, &
 !                         INSERT_VALUES, grid%iphas_loc, ierr)

  call VecGetArrayF90(grid%xx_loc, locpat%xx_loc_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, locpat%porosity_loc_p, ierr)
  call VecGetArrayF90(grid%tor_loc, locpat%tor_loc_p, ierr)
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, locpat%perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, locpat%perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, locpat%perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, locpat%volume_p, ierr)

  call VecGetArrayF90(grid%ithrm_loc, locpat%ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, locpat%icap_loc_p, ierr)
  

  call IMSJacobin_patch(locpat%xx_loc_p, Jac_elem, Jac_ind,locpat,grid,ierr)


! Then set values into matrix
do n= 1, locpat%nlmax
      ng = locpat%nL2G(n)
      if(ng>0)then
	  na1=Jac_ind(n,0)
      !if(n==1) print *, n,na1, Jac_elem(n,0,:,:)
	  blkmat11(:,:)=Jac_elem(n,0,:,:)
      !if(n==1) print *, blkmat11(:,:), Jac_elem(n,0,:,:)
	   call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
       do ind=1,6
         na2=Jac_ind(n,ind)
	     if(na2>=0)then
         blkmat11(:,:)=Jac_elem(n,ind,:,:)  
        ! if(n==1) print *, blkmat11(:,:)   
         call MatSetValuesBlocked(A,1,na1,1,na2,blkmat11,ADD_VALUES,ierr)
         endif
       enddo  
     endif 
enddo  
 


  call VecRestoreArrayF90(grid%xx_loc, locpat%xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, locpat%porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%tor_loc, locpat%tor_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, locpat%perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, locpat%perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, locpat%perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, locpat%volume_p, ierr)

   
    call VecRestoreArrayF90(grid%ithrm_loc, locpat%ithrm_loc_p, ierr)
    call VecRestoreArrayF90(grid%icap_loc, locpat%icap_loc_p, ierr)
!    call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)


 ! call MatCpoy(A,B,ierr)
  
 !call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB, ierr)
   

! call MatView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
! stop

 end subroutine IMSJacobin




 subroutine pflow_IMS_initaccum(grid, locpat)
 
  use translator_ims_module  
  use IMS_patch_module, only : IMSRes_ARCont
  implicit none
  type(pflowGrid) :: grid 
  type(pflow_localpatch_info), intent(inout) :: locpat

  integer :: ierr
  integer*4 :: n,ng
  integer*4 :: i, index_var_begin,index_var_end
  integer*4 :: p1
  integer*4 :: ii1,ii2, iicap

  PetscScalar, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_p(:),&
                          var_p(:), icap_p(:),ithrm_p(:)
  
 !  integer, pointer ::iphase_p(:)
  
  real*8 :: sat_pressure, pvol, satw, tmp, temp  ! Saturation pressure of water.
  real*8 :: dif(1:grid%nphase),res(1:grid%ndof)
 
 
  temp=grid%tref
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
 ! call VecGetArrayF90(grid%var, var_p,ierr)
 
  
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  !print *,'IMS initaccum  Gotten pointers'
 
 
 if (grid%SAMRAI_drive == PETSC_TRUE)then
   ! SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI
 ! var_p =>  locpat%var
 
  !do n = 1, locpat%nlmax
        
   !     iicap=int(icap_p(n))
!		ng =  locpat%nL2G(n)
		!dif(1)= grid%difaq; dif(3)=dif(1)
        !dif(2)= grid%cdiff(int(ithrm_p(n)))

 !  	call pri_var_trans_ims_ninc(yy_p((n-1)*grid%ndof+1:n*grid%ndof), temp,&
 !       grid%scale,grid%nphase, &
 !       iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
 !       grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
 !       grid%pcbetac(iicap),grid%pwrprm(iicap),&
!		var_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),ierr)


 !enddo

  !call VecRestoreArrayF90(grid%var, var_p,ierr)
  !call VecGetArrayF90(grid%var, var_p,ierr)

!---------------------------------------------------------------------------


 ! do n = 1, locpat%nlmax  ! For each local node do...
  !  ng = grid%nL2G(n)   ! corresponding ghost index
 !   p1 = 1 + (n-1)*grid%ndof
!	ng =  locpat%nL2G(n)
!    index_var_begin=(ng-1)*grid%size_var_node+1
!    index_var_end = index_var_begin -1 + grid%size_var_use
!    i = ithrm_p(n)
    
 !    call IMSRes_ARCont(n, var_p(index_var_begin: index_var_end),&
!	  porosity_p(n),volume_p(n),grid%dencpr(i), locpat,grid, Res, 0,ierr)
 

!	accum_p(p1:p1+grid%ndof-1)=Res(:) 

! end do
  
   ! SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI
 else
 var_p =>  locpat%var
 
  do n = 1, locpat%nlmax
        
        iicap=int(icap_p(n))
		ng =  locpat%nL2G(n)
		!dif(1)= grid%difaq; dif(3)=dif(1)
        !dif(2)= grid%cdiff(int(ithrm_p(n)))

   	call pri_var_trans_ims_ninc(yy_p((n-1)*grid%ndof+1:n*grid%ndof), temp,&
        grid%scale,grid%nphase, &
        iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),&
		var_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),ierr)


 enddo

  !call VecRestoreArrayF90(grid%var, var_p,ierr)
  !call VecGetArrayF90(grid%var, var_p,ierr)

!---------------------------------------------------------------------------


  do n = 1, locpat%nlmax  ! For each local node do...
  !  ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
	ng =  locpat%nL2G(n)
    index_var_begin=(ng-1)*grid%size_var_node+1
    index_var_end = index_var_begin -1 + grid%size_var_use
    i = ithrm_p(n)
    
     call IMSRes_ARCont(n, var_p(index_var_begin: index_var_end),&
	  porosity_p(n),volume_p(n),grid%dencpr(i), locpat,grid, Res, 0,ierr)
 

	accum_p(p1:p1+grid%ndof-1)=Res(:) 

 end do
  nullify(var_p)
endif


  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%accum, accum_p, ierr)
!  call VecRestoreArrayF90(grid%var, var_p,ierr)
 
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%icap, icap_p, ierr)

 end subroutine pflow_IMS_initaccum


  subroutine pflow_update_IMS(grid, locpat)
    use translator_ims_module  
    ! use water_eos_module
    implicit none
    type(pflowGrid) :: grid 
    type(pflow_localpatch_info), intent(inout) :: locpat
                     
    integer :: n, ichange,n0, ng
    integer :: ierr,iicap
	PetscScalar, pointer :: xx_p(:),icap_p(:),ithrm_p(:), var_p(:)
    real*8  tmp, temp           

      
  ! if (grid%rk > 0.d0) call Rock_Change(grid)
  ! if(grid%dt<1.D-1*grid%tconv) &
  ! call  translator_IMS_Switching(grid%xx,grid%tref,grid,1,ichange, ierr)
  !print *,'IMS_Update done'
 
   ! if(ichange ==1)then
   temp=grid%tref
   call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
   call VecGetArrayF90(grid%icap,icap_p,ierr)
   call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)  
!   call VecGetArrayF90(grid%var,var_p,ierr)
  


if (grid%SAMRAI_drive == PETSC_TRUE)then
  print *,"SAMRAI SAMRAI SAMRAI SAMRAI"
   ! SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI
 !   do ipatch =1, npatch
 !    do n = 1, grid%locpat(ipatch)%nlmax
 !   iicap = icap_p(n)
 !   n0=(n-1)*grid%ndof
!	ng=locpat%nL2G(n)
!    if(xx_p(n0+3)<0.D0) xx_p(n0+3)=0.D0

!     call pri_var_trans_ims_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),&
!        temp, grid%scale,grid%nphase,&
!        iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
!        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
!        grid%pcbetac(iicap),grid%pwrprm(iicap),&
!		var_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),&
	!	ierr)
  ! enddo
 ! enddo
 
   ! SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI SAMRAI
 else
   var_p=>locpat%var
   do n = 1, locpat%nlmax
    iicap = icap_p(n)
    n0=(n-1)*grid%ndof
	ng=locpat%nL2G(n)
!    if(xx_p(n0+3)<0.D0) xx_p(n0+3)=0.D0

     call pri_var_trans_ims_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),&
        temp, grid%scale,grid%nphase,&
        iicap, grid%sir(1:grid%nphase,iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),&
		var_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use),&
		ierr)
      ! print *,"Var::",var_p((ng-1)*grid%size_var_node+1:(ng-1)*grid%size_var_node+grid%size_var_use)
   enddo
    nullify(var_p) 
 endif
   
   call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
   call VecRestoreArrayF90(grid%icap,icap_p,ierr)
   call VecRestoreArrayF90(grid%ithrm,ithrm_p,ierr)  
   
 !  call VecRestoreArrayF90(grid%var,var_p,ierr)
  
   
  
   if(grid%nphase>1) call pflow_IMS_massbal(grid, locpat)
 ! endif 
!   print *, 'ims_update: end massbal' 

  call VecCopy(grid%xx, grid%yy, ierr)   
   
  call  pflow_ims_initaccum(grid, locpat)
  !print *,'pflow_IMS_initaccum done'
 
 ! call translator_get_output(grid)
 
 
 ! print *,'translator_get_output done'
  ! the output variables should be put into grid%pressure, temp,xmol,sat...
  ! otherwise need to rewrite the pflow_output


 end subroutine pflow_update_ims





end module IMS_module
