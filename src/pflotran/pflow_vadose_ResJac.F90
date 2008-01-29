! introduced grid variables: e_total :: 1 dof
!Vadose translator_module                           c_total :: grid%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (grid%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

               
module Vadose_module

  use pflow_gridtype_module
 ! use pflow_var_module

  private 
  
#include "definitions.h"
!#include "include/petscf90.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
!#ifdef USE_PETSC216
!#include "include/finclude/petscsles.h"
!#endif
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"


! Cutoff parameters
  PetscReal, parameter :: formeps   = 5.D-5
  PetscReal, parameter :: eps       = 1.D-5
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: dfac = 1.D-8
  

  PetscInt,save :: size_var_use 
  PetscInt,save :: size_var_node
  PetscReal, allocatable,save :: Resold_AR(:,:), Resold_FL(:,:)
  PetscReal, pointer, save :: yybc(:,:), vel_bc(:,:)
! Contributions to residual from accumlation/source/Reaction, flux(include diffusion)
  
  
  
   

  public VadoseResidual, VadoseJacobian, pflow_Vadose_initaccum, &
         pflow_update_Vadose,pflow_Vadose_initadj, pflow_Vadose_timecut,&
         pflow_Vadose_setupini, Vadose_Update, Vadose_Update_Reason, &
         pflow_Vadose_bcadj



contains


subroutine pflow_Vadose_timecut(grid)
 
  implicit none
  type(pflowGrid), intent(inout) :: grid
  
 
  PetscReal, pointer :: xx_p(:),yy_p(:)!,var_p(:),iphase_p(:)
  PetscInt :: n,n0,re
  PetscErrorCode :: ierr
  !PetscReal, pointer :: sat(:),xmol(:)

  call VecGetArrayF90(grid%xx, xx_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr)
 ! call VecGetArrayF90(grid%var, var_p, ierr); 
 ! call VecGetArrayF90(grid%iphas, iphase_p, ierr); 

  do n=1, grid%nlmax
    n0=(n-1)*grid%ndof
    do re = 1, grid%ndof
   !xx_p(n0+re)= 0.5D0 * xx_p(n0+re) +.5D0 *yy_p(n0+re)
      xx_p(n0+re)= yy_p(n0+re)
    enddo
  enddo 
  call VecRestoreArrayF90(grid%xx, xx_p, ierr) 
  call VecRestoreArrayF90(grid%yy, yy_p, ierr)
  
  !call VecCopy(grid%xx,grid%yy,ierr)
  !call pflow_Vadose_initaccum(grid)
 
end subroutine pflow_Vadose_timecut
  

subroutine pflow_Vadose_setupini(grid)
  implicit none
  type(pflowGrid), intent(inout) :: grid
  
  PetscReal, pointer :: xx_p(:), iphase_p(:)
  PetscInt :: iln,na,nx,ny,nz,ir
  PetscErrorCode :: ierr
  
  size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
  size_var_node = (grid%ndof + 1) * size_var_use
  
  allocate(Resold_AR(grid%nlmax,grid%ndof))
  allocate(Resold_FL(grid%nconn,grid%ndof))
  allocate(grid%delx(grid%ndof,grid%ngmax))
  grid%delx=0.D0
   
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)
  
  do iln=1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(iln)) <= 0) cycle
    endif

    na = grid%nL2A(iln)
    
!   nz = int(na/grid%nxy) + 1
!   ny = int(mod(na,grid%nxy)/grid%nx) + 1
!   nx = mod(mod(na,grid%nxy),grid%nx) + 1
    
    !compute i,j,k indices from na: note-na starts at 0
    nz = na/grid%nxy + 1
    ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
    nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
    
!   print *,'pflow_Vadose_resjac: ',na,nx,ny,nz
    
    do ir = 1,grid%iregini
      if ((nz>=grid%k1ini(ir)) .and. (nz<=grid%k2ini(ir)) .and.&
          (ny>=grid%j1ini(ir)) .and. (ny<=grid%j2ini(ir)) .and.&
          (nx>= grid%i1ini(ir)) .and. (nx<=grid%i2ini(ir))) then
        iphase_p(iln)=grid%iphas_ini(ir)
        xx_p(1+(iln-1)*grid%ndof:iln*grid%ndof) = grid%xx_ini(:,ir)
            !exit
      endif
    enddo 
  enddo
              
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)

end  subroutine pflow_Vadose_setupini
  

subroutine Vadose_Update_Reason(reason,grid)
  
  implicit none
 
  PetscInt, intent(out) :: reason
  type(pflowGrid), intent(inout) :: grid
  PetscReal, pointer :: xx_p(:),var_p(:),iphase_p(:), yy_p(:) !,r_p(:)
  PetscInt :: n,n0,re
  PetscInt :: re0, iipha
  PetscErrorCode :: ierr
! PetscInt :: index
! PetscReal, pointer :: sat(:),xmol(:)
! PetscReal :: rmax(grid%ndof)

  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  re = 1
 ! call SNESComputeFunction(grid%snes,grid%xx,grid%r,ierr)
 ! do n=1,grid%ndof
 !  call VecStrideNorm(grid%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if (re>0) then
    call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(grid%yy, yy_p, ierr)
    call VecGetArrayF90(grid%var, var_p, ierr); 
    call VecGetArrayF90(grid%iphas, iphase_p, ierr); 
  
    do n = 1,grid%nlmax

      !geh - Ignore inactive cells with inactive materials
      if (associated(grid%imat)) then
        if (grid%imat(grid%nL2G(n)) <= 0) cycle
      endif

      n0=(n-1)* grid%ndof
      !index=(n-1)*size_var_node
      !sat=>var_p(index+2+1:index+2+grid%nphase)
      !den=>var_p(index+2+grid%nphase+1:index+2+2*grid%nphase)
    !xmol=>var_p(index+2+7*grid%nphase+1:index+2+7*grid%nphase + grid%nphase*grid%nspec)    
      iipha=int(iphase_p(n))
     !if(n==3583 .or. n==3587)
   !print *, 'update reson', grid%nlmax, n, iipha, xx_p(n0+1:n0+3)
   !if(xmol(4)>1.0) re=0; goto 1
   !if(xmol(4)<.0) re=0; goto 1
   !if(sat(2) < .0) re=0;goto 1
   ! if(sat(2) > 1.) re=0;goto 1
      select case(iipha)
        case (1)
          if (xx_p(n0 + 3) > 1.0D0) then
            re = 0
            exit
!    goto 111
          endif
          if (xx_p(n0 + 3) < 0.D0) then
            re = 0
            exit
!    goto 111
          endif
     !if(xx_p(n0 + 3) > 1.0D0) xx_p(n0 + 3)=1.D0
     !if(xx_p(n0 + 3) < .0D0) xx_p(n0 + 3)=0.D0
        case (2)
          if (xx_p(n0 + 3) > 1.0D0) then
            re = 0
            exit
!    goto 111
          endif
          if (xx_p(n0 + 3) < 0.D0) then
            re = 0
            exit
!    goto 111
          endif
        case (3)
          if (xx_p(n0 + 3) > 1.D0) then
            re = 0
            exit
!     goto 111
          endif
          if (xx_p(n0 + 3) < 0.D0) then
            re = 0
            exit
!     goto 111
          endif
     !if(xx_p(n0 + 3) > 1.0D0) xx_p(n0 + 3)=1.D0
     !if(xx_p(n0 + 3) < .0D0) xx_p(n0 + 3)=0.D0
      end select  
    enddo
  
!  do n = 1,grid%nlmax
!     n0=(n-1)* grid%ndof
!      
!   if(dabs(xx_p(n0+1)-yy_p(n0+1))>1D6) then
!      re=0;exit
!   endif
   
!   if(dabs(xx_p(n0+2)-yy_p(n0+2))>1D1) then
!      re=0;exit
!   endif
  
  
!   enddo
   ! print *, 'update reason: ',grid%myrank,grid%nlmax,n,re

  !   call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
   !print *,' update reason ba MPI', ierr
    if (re<=0) print *,'Sat or Con out of Region at: ',n,iipha,xx_p(n0+1:n0+3)
    call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(grid%yy, yy_p, ierr)
    call VecRestoreArrayF90(grid%var, var_p, ierr) 
    call VecRestoreArrayF90(grid%iphas, iphase_p, ierr) 
  endif
 ! print *,' update reason', grid%myrank, re,n,grid%nlmax

  
  if (grid%commsize >1) then
    call MPI_ALLREDUCE(re, re0,ONE_INTEGER, MPI_INTEGER,MPI_SUM, &
                       PETSC_COMM_WORLD,ierr)
  !print *,' update reason re'
    !call MPI_BCAST(re0,ONE_INTEGER, MPI_INTEGER, ZERO_INTEGER,PETSC_COMM_WORLD,ierr)
  !print *,' update reason ca'
    if (re0<grid%commsize) re = 0
  endif
  reason = re
  
  if (reason<=0) print *,'Sat or Con out of Region'
  
end subroutine Vadose_Update_Reason

!=======================================================================================
 


subroutine VadoseRes_ARCont(node_no, var_node,por,vol,rock_dencpr, grid, Res_AR,ireac,ierr)

  implicit none

  PetscInt :: node_no
  PetscInt, optional:: ireac
  PetscErrorCode :: ierr
  type(pflowGrid), intent(in) :: grid
  PetscReal, target:: var_node(1:size_var_use)
  PetscReal Res_AR(1:grid%ndof) 
  PetscReal vol,por,rock_dencpr
     
  PetscReal, pointer :: temp, pre_ref   ! 1 dof
  PetscReal, pointer :: sat(:), density(:), amw(:), h(:), u(:), pc(:), kvr(:)         ! nphase dof
  PetscReal, pointer :: xmol(:), diff(:)            ! nphase*nspec
  
  PetscInt :: ibase, m,np, iireac=1
  PetscReal pvol,mol(grid%nspec),eng
  
  if (present(ireac)) iireac=ireac
  pvol=vol*por
  
  ibase=1;                 temp=>var_node(ibase)
  ibase=ibase+1;           pre_ref=>var_node(ibase)
  ibase=ibase+1;           sat=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; amw=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; h=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; u=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr=>var_node(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; xmol=>var_node(ibase:ibase+grid%nphase*grid%nspec-1)
  ibase=ibase+grid%nphase*grid%nspec; diff=>var_node(ibase:ibase+grid%nphase*grid%nspec-1)

  !sumation of component
  mol=0.D0; eng=0.D0
  do np = 1, grid%nphase     
    do m=1, grid%nspec  
      mol(m) = mol(m) + sat(np)*density(np)*xmol(m + (np-1)*grid%nspec)
    enddo
  eng = eng + density(np)*u(np) *sat(np)
  enddo

  mol = mol * pvol
  eng = eng * pvol + (1.D0 - por)* vol * rock_dencpr * temp 
  
! Reaction terms here
  if (iireac>0) then
!H2O
    mol(1)= mol(1) - grid%dt * grid%rtot(node_no,1)
!CO2
    mol(2)= mol(2) - grid%dt * grid%rtot(node_no,2)
!should include related energy change here
  endif
  Res_AR(1:grid%ndof-1)=mol(:)
  Res_AR(grid%ndof)=eng
  nullify(temp, pre_ref, sat, density, amw, h,u, pc,kvr,xmol,diff)       

end subroutine  VadoseRes_ARCont


subroutine VadoseRes_FLCont(nconn_no,area,var_node1,por1,tor1,sir1,dd1,perm1, &
                            Dk1,var_node2,por2,tor2,sir2,dd2,perm2,Dk2,grid, &
                            vv_darcy,Res_FL)

  implicit none
  
  PetscInt :: nconn_no
  type(pflowGrid), intent(inout) :: grid
  PetscReal sir1(1:grid%nphase),sir2(1:grid%nphase)
  PetscReal, target:: var_node1(1:2+7*grid%nphase+2*grid%nphase*grid%nspec)
  PetscReal, target:: var_node2(1:2+7*grid%nphase+2*grid%nphase*grid%nspec)
  PetscReal por1,por2,tor1,tor2,perm1,perm2,Dk1,Dk2,dd1,dd2
  PetscReal vv_darcy(grid%nphase),area
  PetscReal Res_FL(1:grid%ndof) 
     
  PetscReal, pointer :: temp1, pre_ref1   ! 1 dof
  PetscReal, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  PetscReal, pointer :: xmol1(:), diff1(:)            ! 
  
  PetscReal, pointer :: temp2, pre_ref2   ! 1 dof
  PetscReal, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  PetscReal, pointer :: xmol2(:), diff2(:)    
  
  PetscInt :: ibase, m,np, ind
  PetscReal  fluxm(grid%nspec),fluxe, v_darcy,q
  PetscReal uh,uxmol(1:grid%nspec), ukvr,difff,diffdp, DK,Dq
  PetscReal upweight,density_ave,cond, gravity, dphi
  
!  m1=grid%nd1(nc); n1 = grid%nG2L(m1) ! = zero for ghost nodes 
!  print *,'in FLcont'
  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
               
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
               
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+grid%nphase-1)
                           sat2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density1=>var_node1(ibase:ibase+grid%nphase-1)
                           density2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; amw1=>var_node1(ibase:ibase+grid%nphase-1)
                           amw2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; h1=>var_node1(ibase:ibase+grid%nphase-1)
                           h2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; u1=>var_node1(ibase:ibase+grid%nphase-1)
                           u2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc1=>var_node1(ibase:ibase+grid%nphase-1)
                           pc2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr1=>var_node1(ibase:ibase+grid%nphase-1)
                           kvr2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; 
                         xmol1=>var_node1(ibase:ibase+grid%nphase*grid%nspec-1)
                         xmol2=>var_node2(ibase:ibase+grid%nphase*grid%nspec-1)
  ibase=ibase+grid%nphase*grid%nspec;
                         diff1=>var_node1(ibase:ibase+grid%nphase*grid%nspec-1)
                         diff2=>var_node2(ibase:ibase+grid%nphase*grid%nspec-1)

  !print *,' FLcont got pointers' ,var_node1,var_node2,sir1,sir2
  !print *,' tmp=',temp1,temp2
  !print *,'diff=',diff1,diff2
   
  Dq = (perm1 * perm2)/(dd1*perm2 + dd2*perm1)
  diffdp = (por1 *tor1 * por2*tor2) / (dd2*por1*tor1 + dd1*por2*tor2)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy = 0.D0  
  
  do np=1, grid%nphase

! Flow term
    if ((sat1(np) > sir1(np)) .or. (sat2(np) > sir2(np))) then
    
      upweight=dd1/(dd1+dd2)
      if (sat1(np) <eps) then 
        upweight=0.d0
      else if (sat2(np) <eps) then 
        upweight=1.d0
      endif
      density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
    
      gravity = (upweight*density1(np)*amw1(np) + &
                (1.D0-upweight)*density2(np)*amw2(np)) &
                * grid%gravity * grid%delz(nconn_no) * grid%grav_ang(nconn_no)
        
      dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
!    print *,'FLcont  dp',dphi
  ! note uxmol only contains one phase xmol
      if (dphi>=0.D0) then
        ukvr = kvr1(np)
        uh = h1(np)
        uxmol(1:grid%nspec) = xmol1((np-1)*grid%nspec+1:np*grid%nspec)
      else
        ukvr = kvr2(np)
        uh = h2(np)
        uxmol(1:grid%nspec) = xmol2((np-1)*grid%nspec+1:np*grid%nspec)
      endif      
     
   ! print *,'FLcont  uxmol',uxmol
      if (ukvr>floweps) then
        v_darcy= Dq * ukvr * dphi
         !grid%vvl_loc(nconn_no) = v_darcy
        vv_darcy(np) = v_darcy
     
        q = v_darcy * area
          
        do m=1, grid%nspec 
          fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
        enddo

        fluxe = fluxe + q*density_ave*uh 
      endif
    endif 
 !  print *,' FLcont end flow',np
! Diffusion term   
! Note : average rule may not be correct  
    if ((sat1(np) > eps) .and. (sat2(np) > eps)) then
     
      difff = diffdp * 0.25D0*(sat1(np)+sat2(np))*(density1(np)+density2(np))
      do m=1, grid%nspec
        ind=m+(np-1)*grid%nspec
        fluxm(m) = fluxm(m) + difff * .5D0 * &
                   (diff1(ind) + diff2(ind))*(xmol1(ind)-xmol2(ind))
      enddo  
    endif 
  !   print *,' FLcont',np,fluxm,fluxe
  enddo
   
! conduction term
        
  Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)
  cond = Dk*area*(temp1-temp2) 
  fluxe=fluxe + cond
   !      print *,' FLcont heat cond', Dk, cond
  Res_FL(1:grid%ndof-1) = fluxm(:) * grid%dt
  Res_FL(grid%ndof) = fluxe * grid%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 !                           2 R = R - Res_FL  
 !print *,'end FLcont'
 
  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)       
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)       

end subroutine VadoseRes_FLCont

subroutine VadoseRes_FLBCCont(nbc_no,area,var_node1,var_node2,por2,tor2,sir2, &
                              dd1,perm2,Dk2,grid,vv_darcy,Res_FL)
 ! Notice : index 1 stands for BC node
  implicit none
  
  PetscInt :: nbc_no
  type(pflowGrid), intent(inout) :: grid
  PetscReal dd1, sir2(1:grid%nphase)
  PetscReal, target:: var_node1(1:2+7*grid%nphase+2*grid%nphase*grid%nspec)
  PetscReal, target:: var_node2(1:2+7*grid%nphase+2*grid%nphase*grid%nspec)
  PetscReal por2,perm2,Dk2,tor2
  PetscReal vv_darcy(grid%nphase), area
  PetscReal Res_FL(1:grid%ndof) 
     
  PetscReal, pointer :: temp1, pre_ref1   ! 1 dof
  PetscReal, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  PetscReal, pointer :: xmol1(:), diff1(:)            ! 
  
  PetscReal, pointer :: temp2, pre_ref2   ! 1 dof
  PetscReal, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  PetscReal, pointer :: xmol2(:), diff2(:)    
  
  PetscInt :: ibase, m,np, ind, ibc,j
  PetscReal  fluxm(grid%nspec),fluxe, v_darcy,q
  PetscReal uh,uxmol(1:grid%nspec), ukvr,diff,diffdp, DK,Dq
  PetscReal upweight,density_ave,cond,gravity, dphi

  
  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+grid%nphase-1)
                           sat2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; density1=>var_node1(ibase:ibase+grid%nphase-1)
                           density2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; amw1=>var_node1(ibase:ibase+grid%nphase-1)
                           amw2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; h1=>var_node1(ibase:ibase+grid%nphase-1)
                           h2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; u1=>var_node1(ibase:ibase+grid%nphase-1)
                           u2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; pc1=>var_node1(ibase:ibase+grid%nphase-1)
                           pc2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; kvr1=>var_node1(ibase:ibase+grid%nphase-1)
                           kvr2=>var_node2(ibase:ibase+grid%nphase-1)
  ibase=ibase+grid%nphase; 
                         xmol1=>var_node1(ibase:ibase+grid%nphase*grid%nspec-1)
                         xmol2=>var_node2(ibase:ibase+grid%nphase*grid%nspec-1)
  ibase=ibase+grid%nphase*grid%nspec;
                         diff1=>var_node1(ibase:ibase+grid%nphase*grid%nspec-1)    
                         diff2=>var_node2(ibase:ibase+grid%nphase*grid%nspec-1)


  ibc = grid%ibconn(nbc_no)
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy = 0.D0 
   
  select case (grid%ibndtyp(ibc))
    case(1) 
      Dq = perm2 / dd1
      diffdp = por2*tor2/dd1*area
        ! Flow term
      do np=1,grid%nphase
        if ((sat1(np) > sir2(np)) .or. (sat2(np) > sir2(np))) then
          upweight=.5D0
        if (sat1(np) <eps) then 
          upweight=0.d0
        else if (sat2(np) <eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
   
        gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * grid%gravity * grid%delzbc(nbc_no)
       
        dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
   
   
        if (dphi>=0.D0) then
          ukvr = kvr1(np)
          uh = h1(np)
          uxmol(:)=xmol1((np-1)*grid%nspec+1 : np*grid%nspec)
        else
          ukvr = kvr2(np)
          uh = h2(np)
          uxmol(:)=xmol2((np-1)*grid%nspec+1 : np*grid%nspec)
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        !grid%vvl_loc(nbc_no) = v_darcy
          vv_darcy(np) = v_darcy
     
          q = v_darcy * area
          
          do m=1, grid%nspec 
            fluxm(m) = fluxm(m) + q*density_ave*uxmol(m)
          enddo 

          fluxe = fluxe + q*density_ave*uh 
        endif
      endif 
! Diffusion term   
! Note : average rule may not be correct  
      if ((sat1(np) > eps) .and. (sat2(np) > eps)) then
     
        diff = diffdp * 0.25D0*(sat1(np)+sat2(np))*(density1(np)+density2(np))
        do m = 1, grid%nspec
          ind=m+(np-1)*grid%nspec
          fluxm(m) = fluxm(m) + diff * diff2(ind)*( xmol1(ind)-xmol2(ind))
        enddo  
      endif
    enddo
! conduction term
        
    Dk =  Dk2 / dd1
    cond = Dk*area*(temp1-temp2) 
    fluxe=fluxe + cond
 
    Res_FL(1:grid%nspec)=fluxm(:)* grid%dt
    Res_FL(grid%ndof)=fluxe * grid%dt

  case(2)
    if ((dabs(grid%velocitybc(1,nbc_no))+ &
         dabs(grid%velocitybc(2,nbc_no)))>floweps) then
!geh      print *, 'FlowBC :', nbc_no,grid%velocitybc(1,nbc_no), &
!geh               grid%velocitybc(2,nbc_no)

           ! select case(grid%iphasebc(nbc_no))
           ! case(1)
           np = 2
           upweight=.5
            gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * grid%gravity * grid%delzbc(nbc_no)
       
            dphi = pre_ref1 -  pre_ref2  + gravity  
            Dq = perm2 / dd1
             grid%velocitybc(1,nbc_no) = vel_bc(1,nbc_no)
             grid%velocitybc(2,nbc_no) = Dq * kvr2(2) * dphi
          !print *,'FLBC:: ',grid%velocitybc(:,nbc_no), pre_ref1, pre_ref2, temp1, temp2,h1,h2 
           ! case(2)
           !  grid%velocitybc(2,nbc_no) = vel_bc(2,nbc_no)
           !  grid%velocitybc(1,nbc_no) =0.D0
  
           ! case(3)
           !    grid%velocitybc(1,nbc_no) = vel_bc(1,nbc_no)
               !grid%velocitybc(2,nbc_no) = grid%velocitybc(1,nbc_no) *kvr1(2)/kvr1(1) 
            !  grid%velocitybc(2,nbc_no) = grid%velocitybc(1,nbc_no) *kvr1(2)/kvr1(1)
            !  print *, grid%velocitybc(:,nbc_no)!, kvr2(:)
           ! end select     

      do j=1,grid%nphase
        v_darcy = grid%velocitybc(j,nbc_no)
        vv_darcy(j) = grid%velocitybc(j,nbc_no)
!      grid%vvbc(j+(nc-2)*grid%nphase)= grid%velocitybc(j,nc)
      ! note different from 2 phase version

        if (v_darcy >0.d0) then 
          q = v_darcy * density1(j) * area
       !   print *, density1(j), v_darcy
             !q = 0.d0
             !flux = flux - q
          fluxe = fluxe + q  * h1(j) 
          do m=1, grid%nspec
            fluxm(m) = fluxm(m) + q * xmol1(m + (j-1)*grid%nspec)
          enddo 
        else 
          q = v_darcy * density2(j) * area   
          fluxe = fluxe + q  * h2(j) 
          do m=1, grid%nspec
            fluxm(m) = fluxm(m) + q * xmol2(m + (j-1)*grid%nspec)
          enddo 
        endif 
   
      enddo
    endif
    
    Dk =  Dk2 / dd1
    cond = Dk*area*(temp1-temp2) 
    fluxe=fluxe + cond

    Res_FL(1:grid%nspec) = fluxm(:)*grid%dt
    Res_FL(grid%ndof) = fluxe * grid%dt

  case(3)
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
    
    gravity = (upweight*density1(np)*amw1(np) + (1.D0-upweight)*density2(np)*amw2(np)) &
              * grid%gravity * grid%delzbc(nbc_no)
        
    dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
    
   
    if(dphi>=0.D0)then
       ukvr=kvr1(np)
    !density_ave =  density1(np)
    uh=h1(np)
     uxmol(:)=xmol1((np-1)*grid%nspec+1 : np*grid%nspec)
    else
      ukvr=kvr2(np)
     ! density_ave =  density2(np)
     uh=h2(np)
     uxmol(:)=xmol2((np-1)*grid%nspec+1 : np*grid%nspec)
     endif      
     
    if(ukvr*Dq>floweps)then
         v_darcy= Dq * ukvr * dphi
         !grid%vvl_loc(nbc_no) = v_darcy
     vv_darcy(np)=v_darcy
     
       q=v_darcy * area
          
     do m=1, grid%nspec 
       fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
       enddo
           fluxe = fluxe + q*density_ave*uh 
       endif
       endif 
    enddo
   
    Res_FL(1:grid%nspec)=fluxm(:)* grid%dt
    Res_FL(grid%ndof)=fluxe * grid%dt

  end select

  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)       
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)       

end  subroutine VadoseRes_FLBCCont 


subroutine VadoseResidual(snes,xx,r,grid,ierr)

  use water_eos_module
  use Gas_Eos_Module
  use translator_vad_module
  
  implicit none
 
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(pflowGrid), intent(inout) :: grid

 
  PetscErrorCode :: ierr
  PetscInt :: n, ng, nc, nr
  PetscInt :: i, i1, i2, jn, jng
  PetscInt :: m, m1, m2, n1, n2, ip1, ip2, p1, p2
! PetscInt :: t1, t2, c1, c2, s1, s2
  PetscInt :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
! PetscInt :: i1_hencoeff, i2_hencoeff
  PetscInt :: ibc  ! Index that specifies a boundary condition block
! PetscInt :: j, jm1, jm2, jmu, mu 
! PetscReal :: term1, term2, term3


  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
!              ddensity_p(:), ddensity_loc_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), var_p(:),var_loc_p(:) 
                          
               
! PetscReal, pointer :: pc_p(:), pc_loc_p(:),kvr_p(:), kvr_loc_p(:)

  PetscReal, pointer :: iphase_loc_p(:),icap_p(:),iphase_p(:),&
                          icap_loc_p(:),ithrm_loc_p(:),ithrm_p(:)

  PetscInt :: iicap,iiphase,index_var_begin,index_var_end,iicap1,iicap2,np

! PetscReal :: eng, cond, den, eengl,eengg, &
  PetscReal :: dd1, dd2, &
!           fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl, &
!           fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q, &
            pvoldt, voldt, accum, pvol
! PetscReal :: v_darcy,hflx, por1, por2,Dphi,density_ave, 
  PetscReal :: dd, f1, f2, ff, perm1, perm2
! PetscReal :: D0
! PetscReal :: Dq,Dk  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
! PetscReal :: sat_pressure  ! Saturation pressure of water.
  PetscReal :: dw_kg,dw_mol,dif(grid%nphase)
  PetscReal :: tsrc1,qsrc1,csrc1,enth_src_h2o,enth_src_co2, hsrc1 !, qqsrc
  PetscReal :: cw !,cw1,cw2, xxlw,xxla,xxgw,xxga
! PetscReal :: upweight
! PetscReal :: ukvr,uhh,uconc
  PetscReal :: tmp, rho
! PetscReal :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp,visc
  PetscReal :: Res(grid%ndof),vv_darcy(grid%nphase)
! PetscViewer :: viewer
  PetscInt :: ichange
 
  grid%vvlbc=0.D0
  grid%vvgbc=0.D0
  grid%vvl_loc=0.D0
  grid%vvg_loc=0.D0

  if (grid%iphch<=3) then
    call Translator_Vadose_Switching(xx,grid,ZERO_INTEGER,ichange)
    grid%iphch=grid%iphch+1
  endif  
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  do n = 1, grid%nlmax
    if (xx_p((n-1)*grid%ndof+3) < 0.D0)xx_p((n-1)*grid%ndof+3) = 1.D-6
    if (xx_p((n-1)*grid%ndof+3) > 1.D0)xx_p((n-1)*grid%ndof+3) = 1.D0-1.D-6
  enddo
  call VecRestoreArrayF90(xx, xx_p, ierr)

  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%iphas, &
                           INSERT_VALUES, grid%iphas_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%iphas, &
                          INSERT_VALUES, grid%iphas_loc, ierr)

  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas_loc, iphase_loc_p, ierr); CHKERRQ(ierr)

! there is potential possiblity that the pertubation of p may change the direction of pflow.
! once that happens, code may crash, namely wrong derive. 
  do ng = 1, grid%ngmax 
    !ng=grid%nL2G(n)
    iiphase=int(iphase_loc_p(ng))
  
    grid%delx(1,ng)=xx_loc_p((ng-1)*grid%ndof+1)*dfac*1e-3
    grid%delx(2,ng)=xx_loc_p((ng-1)*grid%ndof+2)*dfac
  
  select case (iiphase)
    case (1)
      if (xx_loc_p((ng-1)*grid%ndof+3) < 0.8) then
        grid%delx(3,ng) =  dfac *xx_loc_p((ng-1)*grid%ndof+3)
      else
       grid%delx(3,ng) =  -dfac *xx_loc_p((ng-1)*grid%ndof+3) 
      endif
      if (grid%delx(3,ng) <1D-9 .and. grid%delx(3,ng)>=0.D0) grid%delx(3,ng) = 1.D-9
      if (grid%delx(3,ng) >-1D-9 .and. grid%delx(3,ng)<0.D0) grid%delx(3,ng) = -1.D-9
   case(2)  
      if (xx_loc_p((ng-1)*grid%ndof+3) <0.8) then
        grid%delx(3,ng) =  dfac *xx_loc_p((ng-1)*grid%ndof+3) 
      else
        grid%delx(3,ng) =  -dfac *xx_loc_p((ng-1)*grid%ndof+3) 
      endif 
      if (grid%delx(3,ng) <1D-9 .and. grid%delx(3,ng)>=0.D0) grid%delx(3,ng) = 1.D-9
      if (grid%delx(3,ng) >-1D-9.and. grid%delx(3,ng)<0.D0) grid%delx(3,ng) = -1.D-9
    case(3)
      if (xx_loc_p((ng-1)*grid%ndof+3) <=0.9) then
        grid%delx(3,ng) = dfac *xx_loc_p((ng-1)*grid%ndof+3) 
      else
        grid%delx(3,ng) = -dfac *xx_loc_p((ng-1)*grid%ndof+3) 
      endif 
      
      if (grid%delx(3,ng) <1D-9 .and. grid%delx(3,ng)>=0.D0) grid%delx(3,ng) = 1.D-9
      if (grid%delx(3,ng) >-1D-9 .and. grid%delx(3,ng)<0.D0) grid%delx(3,ng) = -1.D-9
      
      if ((grid%delx(3,ng)+xx_loc_p((ng-1)*grid%ndof+3))>1.D0) then
        grid%delx(3,ng) = (1.D0-xx_loc_p((ng-1)*grid%ndof+3))/1.D5
      endif
      if ((grid%delx(3,ng)+xx_loc_p((ng-1)*grid%ndof+3))<0.D0) then
        grid%delx(3,ng) = xx_loc_p((ng-1)*grid%ndof+3)/1D5
      endif
    end select
  enddo
  
  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)  
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  call VecGetArrayF90(grid%var,var_p,ierr)
  
  
 ! call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)
!------------------------------------------------------ 





!-----  phase properities ---- last time step---
  do n = 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    jn = 1 + (n-1)*grid%nphase
    ng = grid%nL2G(n)
    ii1 = jn  !1+(n-1)*grid%nphase
    ii2 = n*grid%nphase
    iicap = icap_p(n)
    iiphase = iphase_p(n)
    !*****************
    dif(1)= grid%difaq
    dif(2)= grid%cdiff(int(ithrm_p(n)))
  
  !*******************************************
    call pri_var_trans_vad_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),iiphase, &
                                grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use), &
                                grid%itable,ierr,grid%xxphi_co2(n), &
                                grid%dden_co2(n))



    if (grid%ideriv .eq. 1) then
      call pri_var_trans_vad_winc(xx_p((n-1)*grid%ndof+1:n*grid%ndof), &
                                  grid%delx(1:grid%ndof,ng),iiphase, &
                                  grid%scale,grid%nphase,grid%nspec, iicap, dif,&
                                  var_p((n-1)*size_var_node+size_var_use+1:n* &
                                    size_var_node), &
                                  grid%itable,ierr)
    endif
  
!  print *,'var_p',n,iiphase, var_p((n-1)*size_var_node+1:n*size_var_node)              
!   if(n < 5) print *,'pflow_2ph: ',n,grid%ideriv,grid%xxphi_co2(n)
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%iphas,iphase_p,ierr)
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)
  call VecRestoreArrayF90(grid%ithrm,ithrm_p,ierr)
  call VecRestoreArrayF90(grid%var,var_p,ierr)
  ! call VecRestoreArrayF90(grid%iphase,iphase_p,ierr)
  

  call DAGlobalToLocalBegin(grid%da_var_dof,grid%var, &
                            INSERT_VALUES,grid%var_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_var_dof,grid%var,INSERT_VALUES, &
                          grid%var_loc, ierr)
   
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


! End distribute data 
! now assign access pointer to local variables
  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
! call VecGetArrayF90(grid%yy, yy_p, ierr)
 

  ! notice:: here we assume porosity is constant
 
  call VecGetArrayF90(grid%var_loc,var_loc_p,ierr)
  call VecGetArrayF90(grid%yy,yy_p,ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(grid%vl, vl_p, ierr)
  call VecGetArrayF90(grid%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'


  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
  endif

  Resold_AR=0.D0; ResOld_FL=0.D0

!--------------------------------------------------------------------------
! Calculate accumulation term for interior and exterior nodes.
!--------------------------------------------------------------------------
! print *,grid%rtot
  
!  print *, 'Residual  (init):'
!  print *, r_p
  
  r_p = - accum_p

!  print *, 'Residual  (after accum_p):'
!  print *, r_p

  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
    index_var_begin=(ng-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    
    pvol = volume_p(n)*porosity_loc_p(ng)
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt
    iiphase = iphase_loc_p(ng)
    i = ithrm_loc_p(ng)

    accum = 0.d0
    call VadoseRes_ARCont(n, var_loc_p(index_var_begin: index_var_end),&
    porosity_loc_p(ng),volume_p(n),grid%dencpr(i), grid, Res, ONE_INTEGER,ierr)
   
    r_p(p1:p1+grid%ndof-1) = r_p(p1:p1+grid%ndof-1) + Res(1:grid%ndof)
    Resold_AR(n,1:grid%ndof)= Res(1:grid%ndof) 
  enddo

!  print *, 'Residual  (after accum):'
!  print *, r_p

!************************************************************************
! add source/sink terms
 
  do nr = 1, grid%nblksrc
      
    kk1 = grid%k1src(nr) - grid%nzs
    kk2 = grid%k2src(nr) - grid%nzs
    jj1 = grid%j1src(nr) - grid%nys
    jj2 = grid%j2src(nr) - grid%nys
    ii1 = grid%i1src(nr) - grid%nxs
    ii2 = grid%i2src(nr) - grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1 = grid%qsrc(i,nr)
        csrc1 = grid%csrc(i,nr)
        hsrc1 = grid%hsrc(i,nr)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1 = f1*grid%qsrc(i,nr) + f2*grid%qsrc(i-1,nr)
        csrc1 = f1*grid%csrc(i,nr) + f2*grid%csrc(i-1,nr)
        hsrc1 = f1*grid%hsrc(i,nr) + f2*grid%hsrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
   !print *,'pflow2ph : ', grid%myrank,i,grid%timesrc(i,nr), &
   !grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / grid%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / grid%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.

  if(dabs(hsrc1)>1D-20)then 
       do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
             r_p(n*grid%ndof) = r_p(n*grid%ndof) - hsrc1 * grid%dt   
           enddo
          enddo
       enddo
  endif         

    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)

            call wateos_noderiv(tsrc1,var_loc_p((ng-1)*size_var_node+2), &
                                dw_kg,dw_mol,enth_src_h2o,grid%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
              
            r_p((n-1)*grid%ndof + grid%jh2o) = r_p((n-1)*grid%ndof + grid%jh2o) &
                                               - qsrc1 *grid%dt
            r_p(n*grid%ndof) = r_p(n*grid%ndof) - qsrc1*enth_src_h2o*grid%dt
            Resold_AR(n,grid%jh2o)= Resold_AR(n,grid%jh2o) - qsrc1*grid%dt
            Resold_AR(n,grid%ndof)= Resold_AR(n,grid%ndof) - qsrc1 * &
                                                          enth_src_h2o * grid%dt
      
      
      !           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*grid%fmwh2o, &
!           qsrc1
          enddo
        enddo
      enddo
    endif  
    
    if (csrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*grid%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            call ideal_gaseos_noderiv(var_loc_p((ng-1)*size_var_node+2), &
                 tsrc1,grid%scale,rho,enth_src_co2, tmp)
            enth_src_co2 =enth_src_co2 / grid%fmwa

            r_p((n-1)*grid%ndof + grid%jco2) = r_p((n-1)*grid%ndof + grid%jco2) &
                                               - csrc1*grid%dt
            r_p(n*grid%ndof) = r_p(n*grid%ndof) - csrc1 * enth_src_co2 *grid%dt
            Resold_AR(n,grid%jco2)= Resold_AR(n,grid%jco2) - csrc1*grid%dt
            Resold_AR(n,grid%ndof)= Resold_AR(n,grid%ndof) - csrc1 * &
                                                        enth_src_co2 * grid%dt
       !r_p(s1) = r_p(s1) - csrc1

        !   print *,'pflow2ph_co2: ',grid%myrank,nr,n,ng,tsrc1,rho,grid%fmwco2,csrc1
          enddo
        enddo
      enddo
    endif
  
  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
      
    !  do kk = kk1, kk2
     !   do jj = jj1, jj2
       !   do ii = ii1, ii2
          !    n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
           !   ng = grid%nL2G(n)
           !   p1 = 1+(n-1)*grid%ndof
            !  t1 = p1 + 1
           !   c1 = t1 + 1
          !    qqsrc = qsrc1/ddensity_loc_p(ng)
         !     enth_src = hh_loc_p(ng)
        !      r_p(p1) = r_p(p1) - qsrc1
       !       r_p(t1) = r_p(t1) - qsrc1*enth_src
      !        r_p(c1) = r_p(c1) - qqsrc*CCONC_LOC(ng)
     !     enddo
    !    enddo
   !   enddo
  !  endif
  enddo
  !print *,'finished source/sink term'
  
!  print *, 'Residual  (after source/sink):'
!  print *, r_p

!*********************************************************************


 
! stop
!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------
 
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   

    p1 = 1 + (n1-1)*grid%ndof 
    p2 = 1 + (n2-1)*grid%ndof
   
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)


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
!   if(dabs(perm1-1D-15)>1D-20)print *, 'perm1 error', perm1, ip1, n1,n2
!  if(dabs(perm2-1D-15)>1D-20)print *, 'perm2 error', perm2, ip2, n1,n2
    call VadoseRes_FLCont(nc ,grid%area(nc), &
                          var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m1),tor_loc_p(m1), &
                          grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1, &
                          var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m2),tor_loc_p(m2), &
                          grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,grid, &
                          vv_darcy,Res)
    grid%vvl_loc(nc) = vv_darcy(1)
    grid%vvg_loc(nc) = vv_darcy(2)  
    if (n1 > 0) then               ! If the upstream node is not a ghost node...
      do np =1, grid%nphase 
        vl_p(np+(ip1-1)*grid%nphase+3*grid%nphase*(n1-1)) = vv_darcy(np) 
        ! use for print out of velocity
      enddo
    endif
     
    Resold_FL(nc,1:grid%ndof) = Res(1:grid%ndof) 
    
    if (n1>0) then
      r_p(p1:p1+grid%ndof-1) = r_p(p1:p1+grid%ndof-1) + Res(1:grid%ndof)
    endif
   
    if (n2>0) then
      r_p(p2:p2+grid%ndof-1) = r_p(p2:p2+grid%ndof-1) - Res(1:grid%ndof)
    endif


  enddo
   ! print *,'finished NC' 
 
!  print *, 'Residual  (after flux):'
!  print *, r_p

!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', grid%myrank, grid%sgbc    
 
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    if (ng<=0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif

    p1 = 1 + (m-1) * grid%ndof

    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)

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
         grid%xxbc(:,nc)=xx_loc_p((ng-1)*grid%ndof+1:ng*grid%ndof)
         grid%iphasebc(nc) = int(iphase_loc_p(ng))
        if(dabs(grid%velocitybc(1,nc))>1D-20)then
         
            
          if( grid%velocitybc(1,nc)>0) then
             grid%xxbc(1:2,nc)= yybc(1:2,nc)
          endif     
        endif    
            
            
         
         
      case(3) 
     !  grid%xxbc((nc-1)*grid%ndof+1)=grid%pressurebc(2,ibc)
        grid%xxbc(2:grid%ndof,nc) = xx_loc_p((ng-1)*grid%ndof+2: ng*grid%ndof)
        grid%iphasebc(nc)=int(iphase_loc_p(ng))
    end select

! print *,'2ph bc',grid%myrank,nc,m,ng,ibc,grid%ibndtyp(ibc),grid%pressurebc(:,ibc), &
! grid%tempbc(ibc),grid%sgbc(ibc),grid%concbc(ibc),grid%velocitybc(:,ibc)

!   if(grid%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   grid%ibndphaseRate(ibc) 
   !   grid%ibndconc(ibc)    ! 
   !   grid%tempbc(ibc)      !1 elements 
   !   grid%pressurebc(ibc)  !nphase elements
!      endif
   
    
    iicap=int(icap_loc_p(ng))  
!      print *,'pflow_2pha_bc: ',grid%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      grid%nconnbc,grid%ibndtyp(ibc),grid%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,grid%ideriv,grid%ibndtyp(ibc),grid%density_bc,&
!   grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
       
        !*****************
    dif(1)= grid%difaq
    dif(2)= grid%cdiff(int(ithrm_loc_p(ng)))
    !*******************************************

!    print *, nc
!    print *, 'xxbc: ', grid%xxbc(:,nc)
!    print *, 'iphasebc: ', grid%iphasebc(nc)
!    print *, 'icaptype: ', grid%icaptype(iicap)
!    print *, 'sir: ', grid%sir(1:grid%nphase,iicap)
!    print *, 'lambda: ', grid%lambda(iicap)
!    print *, 'alpha: ', grid%alpha(iicap)
!    print *, 'pckrm: ', grid%pckrm(iicap)
!    print *, 'pwrprm: ', grid%pwrprm(iicap)
!    print *, 'varbc: ', grid%varbc(1:size_var_use)
!    print *, 'xxphi_co2_bc: ', grid%xxphi_co2_bc(nc)
!    print *
  
    call pri_var_trans_vad_ninc(grid%xxbc(:,nc),grid%iphasebc(nc),&
                                grid%scale,grid%nphase,grid%nspec, iicap, dif,&
                                grid%varbc(1:size_var_use),grid%itable,ierr, &
                                grid%xxphi_co2_bc(nc),cw)
   
     !print *,"bc: ", nc, grid%xxbc(:,nc), grid%varbc(1:size_var_use)
    call VadoseRes_FLBCCont(nc,grid%areabc(nc),grid%varbc(1:size_var_use), &
                            var_loc_p((ng-1)*size_var_node+1:(ng-1)* &
                            size_var_node+size_var_use),porosity_loc_p(ng), &
                            tor_loc_p(ng),grid%sir(1:grid%nphase,iicap), &
                            grid%distbc(nc),perm1,D2, grid, vv_darcy,Res)
    grid%vvlbc(nc) = vv_darcy(1)
    grid%vvgbc(nc) = vv_darcy(2) 
    r_p(p1:p1-1+grid%ndof)= r_p(p1:p1-1+grid%ndof) - Res(1:grid%ndof)
    ResOld_AR(m,1:grid%ndof) = ResOld_AR(m,1:grid%ndof) - Res(1:grid%ndof)
   
   
       !print *, ' boundary index', nc,ng,ibc,grid%ibndtyp(ibc)
       !print *,'        xxbc', grid%iphasebc(nc), grid%xxbc(:,nc),res
     !print *, '       var', grid%varbc
  !   print *, ' P  T   C   S  ', grid%pressurebc(1,ibc),grid%tempbc(ibc), &
    !                               grid%concbc(ibc),grid%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),grid%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,grid%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
!print *,grid%density_bc,grid%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc

  enddo
!  print *,'finished BC'

!  print *, 'Residual  (after bc flux):'
!  print *, r_p

  if (grid%use_isoth==PETSC_TRUE) then
    do n = 1, grid%nlmax  ! For each local node do...

      !geh - Ignore inactive cells with inactive materials
      if (associated(grid%imat)) then
        if (grid%imat(grid%nL2G(n)) <= 0) cycle
      endif

      ng = grid%nL2G(n)   ! corresponding ghost index
      p1 = 3 + (n-1)*grid%ndof
      r_p(p1)=xx_loc_p(2 + (ng-1)*grid%ndof)-yy_p(p1-1)
    enddo
  endif

!geh - Zero out rows in matrix and residual entries for inactive cells
  if (associated(grid%imat)) then
    do n = 1, grid%nlmax
      if (grid%imat(grid%nL2G(n)) <= 0) then
        p1 = (n-1)*grid%ndof
        r_p(p1+1:p1+grid%ndof) = 0.d0
      endif
    enddo
  endif

!  print *, 'Residual  (final):'
!  print *, r_p

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr)
  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%accum, accum_p, ierr)
  call VecRestoreArrayF90(grid%var_loc,var_loc_p,ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(grid%vl, vl_p, ierr)
  call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)
  if (grid%rk > 0.d0) then
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif
!#define DEBUG_GEH
#ifdef DEBUG_GEH 
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'residual.out',viewer,ierr)
 call VecView(r,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

  !print *,'XX ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
 !print *,'Residual ::...........'; call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)
 
 !print *,'finished VadoseResidual'
end subroutine VadoseResidual
                
! --------------------------------------------------------------------- 

subroutine VadoseJacobian(snes,xx,A,B,flag,grid,ierr)
       
  use water_eos_module
  use gas_eos_module
  use translator_vad_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(pflowGrid), intent(inout) :: grid
   ! PetscInt, intent(inout) :: flag
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: n, ng, nc,nvar,neq,nr
  PetscInt :: i1, i2, jng, i
  PetscInt :: kk,ii1,jj1,kk1,ii2,jj2,kk2  
  PetscInt :: m, m1, m2, n1, n2, ip1, ip2 
  PetscInt :: p1,p2 !,t1,t2,c1,c2,s1,s2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
! PetscInt :: j, jn, jm1, jm2,jmu, mu
  
! PetscReal :: v_darcy, q
  PetscReal :: dum1, dum2

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), phis_p(:),  tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  PetscInt :: iicap,ii,jj,iiphas,iiphas1,iiphas2,iicap1,iicap2
  PetscInt :: index_var_begin, index_var_end
! PetscInt :: ibc_hencoeff
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho!,dddt,dddp,fg,dfgdp,&
!           dfgdt,eng,dhdt,dhdp,visc,dvdt,dvdp
! PetscReal :: cond, gravity,  acc,  density_ave, 
  PetscReal :: vv_darcy(grid%nphase),&
            voldt, pvoldt
! PetscReal :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
!           flux, fluxh, fluxv, difff, diffg, diffl
  PetscReal :: ff,dif(1:grid%nphase)
  PetscReal :: tsrc1,qsrc1,csrc1
  PetscReal :: dd1, dd2, dd, f1, f2 !, den
! PetscReal :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
! PetscReal :: por1, por2
  PetscReal :: perm1, perm2
! PetscReal :: qu_rate, p_vapor,sat_pressure_t
! PetscReal :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
! PetscReal :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! PetscReal :: sat_pressure  ! Saturation pressure of water.
! PetscReal :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu,sat_ave
  PetscReal :: ra(1:grid%ndof,1:2*grid%ndof)  
! PetscReal :: uhh, uconc, ukvr
  PetscReal :: tmp
! PetscReal :: upweight,m1weight,m2weight,mbweight,mnweight
  PetscReal :: delxbc(1:grid%ndof)
  PetscReal :: blkmat11(1:grid%ndof,1:grid%ndof), &
            blkmat12(1:grid%ndof,1:grid%ndof),&
            blkmat21(1:grid%ndof,1:grid%ndof),&
            blkmat22(1:grid%ndof,1:grid%ndof)
  PetscReal :: ResInc(1:grid%nlmax, 1:grid%ndof, 1:grid%ndof),res(1:grid%ndof)  
  PetscReal :: max_dev  
  PetscInt ::  na1,na2
! PetscViewer :: viewer
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

 ! print *,'*********** In Jacobian ********************** '
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

  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%tor_loc, tor_loc_p, ierr)
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(grid%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(grid%var_loc, var_loc_p, ierr)

 !print *,' In mph Jacobian ::  got pointers '
! ********************************************************************

! Accumulation terms

  ResInc=0.D0
  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ng = grid%nL2G(n)   !get ghosted index
    
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    iiphas=iphase_loc_p(ng)
 ! pressure equation    
    do nvar=1, grid%ndof
   
      index_var_begin=(ng-1)*size_var_node+nvar*size_var_use+1
      index_var_end = index_var_begin -1 + size_var_use

      call VadoseRes_ARCont(n, var_loc_p(index_var_begin : index_var_end), &
                            porosity_loc_p(ng),volume_p(n), &
                            grid%dencpr(int(ithrm_loc_p(ng))),grid, Res,ONE_INTEGER,ierr)
      
      ResInc(n,:,nvar) = ResInc(n,:,nvar) + Res(:)
    enddo
  enddo
! print *,' Mph Jaco Finished accum terms'
! Source / Sink term

  do nr = 1, grid%nblksrc
      
    kk1 = grid%k1src(nr) - grid%nzs
    kk2 = grid%k2src(nr) - grid%nzs
    jj1 = grid%j1src(nr) - grid%nys
    jj2 = grid%j2src(nr) - grid%nys
    ii1 = grid%i1src(nr) - grid%nxs
    ii2 = grid%i2src(nr) - grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1 = grid%qsrc(i,nr)
        csrc1 = grid%csrc(i,nr)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1 = f1*grid%qsrc(i,nr) + f2*grid%qsrc(i-1,nr)
        csrc1 = f1*grid%csrc(i,nr) + f2*grid%csrc(i-1,nr)
        goto 10
      endif
    enddo
10  continue
    
   !print *,'pflow2ph : ', grid%myrank,i,grid%timesrc(i,nr), &
   !grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / grid%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / grid%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
    
    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)
            
            do nvar=1,grid%ndof      
              call wateos_noderiv(tsrc1,var_loc_p((ng-1)*size_var_node+nvar* &
                                                         size_var_use+2), &
                                  dw_kg,dw_mol,enth_src_h2o,grid%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
              ResInc(n,grid%jh2o,nvar) = ResInc(n,grid%jh2o,nvar) - qsrc1 * &
                                                                    grid%dt
              ResInc(n,grid%ndof,nvar) = ResInc(n,grid%ndof,nvar) - qsrc1 * &
                                                        enth_src_h2o * grid%dt

      
      
      !           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*grid%fmwh2o, &
!           qsrc1
            enddo
          enddo
        enddo
      enddo
    endif  
    
    if (csrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*grid%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            do nvar=1,grid%ndof     
              call ideal_gaseos_noderiv(var_loc_p((ng-1)*size_var_node+nvar* &
                                                         size_var_use+2), &
                                        tsrc1,grid%scale,rho,enth_src_co2,tmp)
              enth_src_co2 = enth_src_co2 / grid%fmwa
    
              ResInc(n,grid%jco2,nvar)=  ResInc(n,grid%jco2,nvar) - csrc1 * &
                                                                    grid%dt
              ResInc(n,grid%ndof,nvar)=  ResInc(n,grid%ndof,nvar) - csrc1 * &
                                                           enth_src_co2*grid%dt

          !  Res_AR(n,grid%jco2)= Res_AR(n,grid%jco2) - csrc1
    !  Res_AR(n,grid%ndof)= Res_AR(n,grid%ndof) - csrc1 * enth_src_co2
       !r_p(s1) = r_p(s1) - csrc1

!           print *,'pflow2ph_co2: ',nr,n,ng,tsrc1,rho,grid%fmwco2,csrc1
            enddo
          enddo
        enddo
      enddo
    endif
  enddo  
  
  ! print *,' Mph Jaco Finished source terms'
! Contribution from BC
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    if (ng<=0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif
  
    p1 = 1 + (m-1) * grid%ndof
       
    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)
     
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
        delxbc =0.D0
      case(2)
        ! solve for pb from Darcy's law given qb /= 0
        grid%xxbc(:,nc) = xx_loc_p((ng-1)*grid%ndof+1: ng*grid%ndof)
        grid%iphasebc(nc) = int(iphase_loc_p(ng))
        delxbc = grid%delx(1:grid%ndof,ng)
        
  
       if(dabs(grid%velocitybc(1,nc))>1D-20)then
         if( grid%velocitybc(1,nc)>0) then
             grid%xxbc(1:2,nc)= yybc(1:2,nc)
             delxbc(1:2)=0.D0
          endif     
        endif    
            

         
      case(3) 
        !    grid%xxbc(1,nc)=grid%pressurebc(2,ibc)
        grid%xxbc(2:grid%ndof,nc) = xx_loc_p((ng-1)*grid%ndof+2:ng*grid%ndof)
        grid%iphasebc(nc) = int(iphase_loc_p(ng))
        delxbc(1) = 0.D0
        delxbc(2:grid%ndof) = grid%delx(2:grid%ndof,ng)
    end select

! print *,'2ph bc',grid%myrank,nc,m,ng,ibc,grid%ibndtyp(ibc),grid%pressurebc(:,ibc), &
! grid%tempbc(ibc),grid%sgbc(ibc),grid%concbc(ibc),grid%velocitybc(:,ibc)

!   if(grid%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   grid%ibndphaseRate(ibc) 
   !   grid%ibndconc(ibc)    ! 
   !   grid%tempbc(ibc)      !1 elements 
   !   grid%pressurebc(ibc)  !nphase elements
!      endif
   
    iicap = int(icap_loc_p(ng))     
       
!      print *,'pflow_2pha_bc: ',grid%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      grid%nconnbc,grid%ibndtyp(ibc),grid%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,grid%ideriv,grid%ibndtyp(ibc),grid%density_bc,&
!   grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
        !*****************
    dif(1) = grid%difaq
    dif(2) = grid%cdiff(int(ithrm_loc_p(ng)))
    !*******************************************

  !  print *,' Mph Jaco BC terms: finish setup'
  ! here should pay attention to BC type !!!
    call pri_var_trans_vad_ninc(grid%xxbc(:,nc),grid%iphasebc(nc), &
                                grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                grid%varbc(1:size_var_use),grid%itable,ierr, &
                                dum1, dum2)
  
    call pri_var_trans_vad_winc(grid%xxbc(:,nc),delxbc,grid%iphasebc(nc), &
                                grid%scale,grid%nphase,grid%nspec, iicap,dif(1:grid%nphase),&
                                grid%varbc(size_var_use+1:(grid%ndof+1)* &
                                  size_var_use), &
                                grid%itable,ierr)
            
!    print *,' Mph Jaco BC terms: finish increment'
    do nvar=1,grid%ndof
   
      call VadoseRes_FLBCCont(nc,grid%areabc(nc), &
                              grid%varbc(nvar*size_var_use+1:(nvar+1)* &
                                size_var_use), &
                              var_loc_p((ng-1)*size_var_node+nvar* &
                                size_var_use+1:(ng-1)*size_var_node+nvar* &
                                size_var_use+size_var_use), &
                              porosity_loc_p(ng),tor_loc_p(ng), &
                              grid%sir(1:grid%nphase,iicap), &
                              grid%distbc(nc),perm1,D2,grid,vv_darcy,Res)
    
      ResInc(m,1:grid%ndof,nvar) = ResInc(m,1:grid%ndof,nvar) - Res(1:grid%ndof)
    enddo
 !  print *,' Mph Jaco BC terms: finish comp'
    !   print *, ' boundary index', nc,ng,ibc,grid%ibndtyp(ibc)
    !   print *, ' P  T   C   S  ', grid%pressurebc(1,ibc),grid%tempbc(ibc), &
    !                               grid%concbc(ibc),grid%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),grid%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,grid%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
!print *,grid%density_bc,grid%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc

  enddo
  ! print *,' Mph Jaco Finished BC terms'

  do n= 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ra=0.D0
    ng = grid%nL2G(n)
    na1= grid%nG2N(ng)
   ! Remember, the matrix index starts from (0,0)
    p1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
   
    max_dev = 0.D0
    do neq=1, grid%ndof
      do nvar=1, grid%ndof
        ra(neq,nvar) = ResInc(n,neq,nvar)/grid%delx(nvar,ng) - &
                       ResOld_AR(n,neq)/grid%delx(nvar,ng)
        if (max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
   
      enddo      
    enddo
    if (grid%use_isoth==PETSC_TRUE) then
      ra(:,2)=0.D0
      ra(3,1:grid%ndof)=0.D0
      ra(3,2)=1.D0
    endif     
   
    if (max_dev<1D-5) then
      print *,'Mph Jaco max dev = ', max_dev
    endif
  
    if (grid%iblkfmt == 0) then
      p1=(na1)*grid%ndof
      do ii=0,grid%ndof-1
        do jj=0,grid%ndof-1
          call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
        enddo
      enddo
    else
      blkmat11=ra(1:grid%ndof,1:grid%ndof)
      call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
    endif
         
  enddo
  
!   print *,' Mph Jaco Finished one node terms'
! -----------------------------contribution from transport----------------------

 !print *,'phase cond: ',iphase_loc_p
  ResInc=0.D0
  do nc = 1, grid%nconn  ! For each interior connection...
    ra = 0.D0
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   
    na1 = grid%nG2N(m1)
    na2 = grid%nG2N(m2)
    !print *, grid%myrank,nc,m1,m2,n1,n2,na1,na2
    p1 =  (m1-1)*grid%ndof
    p2 =  (m2-1)*grid%ndof
   
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)
    
    iiphas1 = iphase_loc_p(m1)
    iiphas2 = iphase_loc_p(m2)


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
    iicap1 = int(icap_loc_p(m1))
    iicap2 = int(icap_loc_p(m2))
 
  ! do neq = 1, grid%ndof
    do nvar = 1, grid%ndof
    
      call VadoseRes_FLCont(nc ,grid%area(nc), &
                            var_loc_p((m1-1)*size_var_node+nvar* &
                              size_var_use+1:(m1-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,grid, &
                            vv_darcy,Res)

      ra(:,nvar)= Res(:)/grid%delx(nvar,m1)-ResOld_FL(nc,:)/grid%delx(nvar,m1)
       
      call VadoseRes_FLCont(nc,grid%area(nc), &
                            var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            grid%sir(1:grid%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+nvar* &
                              size_var_use+1:(m2-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            grid%sir(1:grid%nphase,iicap2),dd2,perm2,D2,grid, &
                            vv_darcy,Res)
 
      ra(:,nvar+grid%ndof)= Res(:)/grid%delx(nvar,m2)-ResOld_FL(nc,:)/grid%delx(nvar,m2)
   
    enddo
  
   !   print *,' Mph Jaco Finished NC terms'
  
  ! enddo
   
    if (grid%use_isoth==PETSC_TRUE) then
      ra(3,1:2*grid%ndof)=0.D0
      ra(:,2)=0.D0
      ra(:,2+grid%ndof)=0.D0
    endif   
 
    if (grid%iblkfmt == 1) then
      blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
    endif
    p1=(na1)*grid%ndof;p2=(na2)*grid%ndof
    do ii=0,grid%ndof-1
      do jj=0,grid%ndof-1
        if (n1>0) then
          if (grid%iblkfmt == 0) then
            call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
          else
            blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
          endif
        endif
        if (n2>0) then
          if (grid%iblkfmt == 0) then
            call MatSetValue(A,p2+ii,p1+jj,-ra(ii+1,jj+1),ADD_VALUES,ierr)
          else
            blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
          endif
        endif
      enddo
   
      do jj=grid%ndof,2*grid%ndof-1
        if (n1>0) then
          if (grid%iblkfmt == 0) then
            call MatSetValue(A,p1+ii,p2+jj-grid%ndof,ra(ii+1,jj+1), &
                             ADD_VALUES,ierr)
          else
            blkmat12(ii+1,jj-grid%ndof+1) = blkmat12(ii+1,jj-grid%ndof+1) + &
                                                             ra(ii+1,jj+1)
          endif
        endif
        if (n2>0) then
          if (grid%iblkfmt == 0) then
            call MatSetValue(A,p2+ii,p2+jj-grid%ndof,-ra(ii+1,jj+1), &
                             ADD_VALUES,ierr)
          else
            blkmat22(ii+1,jj-grid%ndof+1) =  blkmat22(ii+1,jj-grid%ndof+1) - &
                                                              ra(ii+1,jj+1)
          endif
        endif
      enddo
    enddo
  
    if (grid%iblkfmt /= 0) then
      if (n1>0) call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
      if (n2>0) call MatSetValuesBlocked(A,1,na2,1,na2,blkmat22,ADD_VALUES,ierr)
      if (n1>0) call MatSetValuesBlocked(A,1,na1,1,na2,blkmat12,ADD_VALUES,ierr)
      if (n2>0) call MatSetValuesBlocked(A,1,na2,1,na1,blkmat21,ADD_VALUES,ierr)
    endif
!print *,'accum r',ra(1:5,1:8)   
 !print *,'devq:',nc,q,dphi,devq(3,:)
  enddo
  ! print *,' Mph Jaco Finished Two node terms'
  
  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(grid%var_loc, var_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)

  if (grid%rk > 0.d0) then
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

!geh - Zero out rows in matrix and residual entries for inactive cells
  if (associated(grid%imat)) then
    do n = 1, grid%nlmax
      ng = grid%nL2G(n)
      if (grid%imat(ng) <= 0) then
        p1=(ng-1)*grid%ndof
        do ii=0,grid%ndof-1
          call MatZeroRowsLocal(A,1,p1+ii,1.d0,ierr)
        enddo
      endif
    enddo
  endif

  !B = A
  !call MatCopy(A,B,ierr)
  
 !call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB, ierr)

#ifdef DEBUG_GEH    
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
! stop

end subroutine VadoseJacobian




subroutine pflow_Vadose_initaccum(grid)
 
  use translator_vad_module  
 
  implicit none
  
  type(pflowGrid) :: grid 

 
  PetscErrorCode :: ierr
  PetscInt :: n
  PetscInt :: i,index_var_begin,index_var_end
  PetscInt :: p1
! PetscInt :: ii1,ii2
  PetscInt :: iicap,iiphase

  PetscReal, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_p(:),&
                          var_p(:), icap_p(:),iphase_p(:),ithrm_p(:)
  
 !PetscInt, pointer ::iphase_p(:)
  
  PetscReal :: pvol, satw ! Saturation pressure of water.
  PetscReal :: dif(1:grid%nphase),res(1:grid%ndof)
  
! PetscReal :: sat_pressure
 
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
  call VecGetArrayF90(grid%var, var_p,ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  !print *,'Vadoseinitaccum  Gotten pointers'
 
  do n = 1, grid%nlmax
        
    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    iicap=int(icap_p(n))
    iiphase = int(iphase_p(n))
    dif(1)= grid%difaq
    dif(2)= grid%cdiff(int(ithrm_p(n)))

    call pri_var_trans_vad_ninc(yy_p((n-1)*grid%ndof+1:n*grid%ndof),iiphase,&
                                grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                size_var_node+size_var_use),grid%itable,ierr, &
                                satw, pvol)
  enddo

  call VecRestoreArrayF90(grid%var, var_p,ierr)
  call VecGetArrayF90(grid%var, var_p,ierr)

!---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

  !  ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
    index_var_begin=(n-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    i = ithrm_p(n)
    
    call VadoseRes_ARCont(n,var_p(index_var_begin:index_var_end), &
                          porosity_p(n),volume_p(n),grid%dencpr(i),grid,Res, &
                          ZERO_INTEGER,ierr)
 

    accum_p(p1:p1+grid%ndof-1)=Res(:) 

   !print *, 'init m accum ', n,  Res 

! print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec),kvr_p(jn),kvr_p(jn+1)

  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%accum, accum_p, ierr)
  call VecRestoreArrayF90(grid%var, var_p,ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%icap, icap_p, ierr)

end subroutine pflow_Vadose_initaccum


subroutine pflow_update_Vadose(grid)

  use translator_vad_module
!  use Condition_module_old
   ! use water_eos_module
  implicit none

  type(pflowGrid) :: grid 
    
! PetscInt :: ichange
  PetscInt :: n,n0
  PetscErrorCode :: ierr
  PetscInt :: iicap,iiphase
  PetscReal, pointer :: xx_p(:),icap_p(:),ithrm_p(:),iphase_p(:), var_p(:)
  PetscReal dif(1:grid%nphase), dum1, dum2           

      
  if (associated(grid%imat)) then
  !  call UpdateBoundaryConditions(grid)
    call pflow_Vadose_bcadj(grid)
  endif

  ! if (grid%rk > 0.d0) call Rock_Change(grid)
  ! call  Translator_vadose_Switching(grid%xx,grid,1,ichange)
  !print *,'Vadose_Update done'
 
   ! if(ichange ==1)then
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)  
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  call VecGetArrayF90(grid%var,var_p,ierr)

  do n = 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    iicap = icap_p(n)
    iiphase = iphase_p(n)
    n0 = (n-1)*grid%ndof
    if (xx_p(n0+3)<0.D0) xx_p(n0+3)=1.D-6

    !*****************
    dif(1) = grid%difaq
    dif(2) = grid%cdiff(int(ithrm_p(n)))
    !*******************************************
    call pri_var_trans_vad_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),iiphase, &
                                grid%scale,grid%nphase,grid%nspec, iicap, dif,&
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use),&
                                grid%itable,ierr, dum1, dum2)

  enddo
   
  call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)
  call VecRestoreArrayF90(grid%ithrm,ithrm_p,ierr)  
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(grid%var,var_p,ierr)
   
  if (grid%nphase>1) call translator_Vadose_massbal(grid)
 ! endif 

  call VecCopy(grid%xx, grid%yy, ierr)   
  call VecCopy(grid%iphas, grid%iphas_old, ierr)   
   
  call  pflow_Vadose_initaccum(grid)
    !print *,'pflow_Vadose_initaccum done'
  call translator_vad_get_output(grid)
 ! print *,'translator_get_output done'
  ! the output variables should be put into grid%pressure, temp,xmol,sat...
  ! otherwise need to rewrite the pflow_output

end subroutine pflow_update_Vadose





subroutine pflow_Vadose_initadj(grid)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  use translator_vad_module  

  implicit none

  type(pflowGrid) :: grid 

 
  PetscErrorCode :: ierr
  PetscInt :: n, nc
  PetscInt :: ibc,jn
  PetscInt :: m
  PetscInt :: ii1,ii2,iicap
  PetscInt :: iiphase,iithrm
 

  PetscReal, pointer :: xx_p(:),var_p(:)
                           

  PetscReal, pointer ::iphase_p(:), ithrm_p(:),icap_p(:)
  
  PetscReal  dif(grid%nphase), dum1, dum2
  
! PetscReal :: temp1
!  PetscReal, parameter :: Rg=8.31415D0

  call VecGetArrayF90(grid%icap, icap_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
  call VecGetArrayF90(grid%xx, xx_p, ierr)
  call VecGetArrayF90(grid%var, var_p, ierr)
! print *,'initadj gotten pointers' 


  do n = 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(grid%imat)) then
      if (grid%imat(grid%nL2G(n)) <= 0) cycle
    endif

    jn = 1 + (n-1)*grid%nphase
    ii1=1+(n-1)*grid%nphase; ii2=n*grid%nphase
    iicap=int(icap_p(n))
        
    iiphase = iphase_p(n)
        !*****************
    dif(1)= grid%difaq
    dif(2)= grid%cdiff(int(ithrm_p(n)))
    !*******************************************
    call pri_var_trans_vad_ninc(xx_p((n-1)*grid%ndof+1:n*grid%ndof),iiphase, &
                                grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use), &
                                grid%itable,ierr, dum1, dum2)
   
    if (translator_check_phase_cond_vad(iiphase, &
                                        var_p((n-1)*size_var_node+1:(n-1)* &
                                        size_var_node+size_var_use), &
                                        grid%nphase,grid%nspec) /= 1 ) then
      print *," Wrong internal node init...  STOP!!!"
      stop    
    endif 
  enddo


  allocate(yybc(grid%ndof,grid%nconnbc))
  allocate(vel_bc(grid%nphase, grid%nconnbc))
  yybc =grid%xxbc
  vel_bc = grid%velocitybc
  
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       
    if (m<0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif

    ibc = grid%ibconn(nc)
       
!      print *,'initadj_bc',nc,ibc,grid%ibndtyp(ibc),grid%nconnbc

    if (grid%ibndtyp(ibc)==1) then
      iicap=int(icap_p(m))
      iithrm=int(ithrm_p(m)) 
      dif(1)= grid%difaq
      dif(2)= grid%cdiff(iithrm)
      call pri_var_trans_vad_ninc(grid%xxbc(:,nc),grid%iphasebc(nc), &
                                  grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                  grid%varbc(1:size_var_use),grid%itable,ierr, &
                                  dum1, dum2)
     ! print *, 'yybc', grid%varbc(1:size_var_use)
      if (translator_check_phase_cond_vad(grid%iphasebc(nc), &
                                          grid%varbc(1:size_var_use), &
                                          grid%nphase,grid%nspec) /=1) then
        print *," Wrong bounday node init...  STOP!!!", grid%xxbc(:,nc)
      
        print *,grid%varbc
        stop    
      endif 
    
               
              
    endif

    if (grid%ibndtyp(ibc)==2) then
  
       yybc(1,nc)= grid%velocitybc(2,nc)
       grid%velocitybc(2,nc) =0.D0
       vel_bc(1,nc) = grid%velocitybc(1,nc)

    endif 

  enddo

  call VecRestoreArrayF90(grid%icap, icap_p, ierr)
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  call VecRestoreArrayF90(grid%var, var_p, ierr)
  !print *,kgjkdf
  
  !need to record initial BC input
  !print *, 'yybc',yybc,grid%varbc
  !call VecCopy(grid%iphas,grid%iphas_old,ierr)
   
end subroutine pflow_Vadose_initadj

!geh We need to be able to update transient boundary conditions.  The
!    below provides the functionality above, but only for bcs.
subroutine pflow_Vadose_bcadj(grid)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  use translator_vad_module  

  implicit none

  type(pflowGrid) :: grid 

 
  PetscErrorCode :: ierr
  PetscInt :: nc, ibc, m, iicap, iithrm                           

  PetscReal, pointer :: ithrm_p(:),icap_p(:)
  
  PetscReal  dif(grid%nphase), dum1, dum2
  
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       
    if (m<0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif

    ibc = grid%ibconn(nc)
       
!      print *,'initadj_bc',nc,ibc,grid%ibndtyp(ibc),grid%nconnbc

    if (grid%ibndtyp(ibc)==1) then
      iicap=int(icap_p(m))
      iithrm=int(ithrm_p(m)) 
      dif(1)= grid%difaq
      dif(2)= grid%cdiff(iithrm)
      call pri_var_trans_vad_ninc(grid%xxbc(:,nc),grid%iphasebc(nc), &
                                  grid%scale,grid%nphase,grid%nspec, iicap, dif, &
                                  grid%varbc(1:size_var_use),grid%itable,ierr, &
                                  dum1, dum2)
      
      if (translator_check_phase_cond_vad(grid%iphasebc(nc), &
                                          grid%varbc(1:size_var_use), &
                                          grid%nphase,grid%nspec) /=1) then
        print *," Wrong bounday node init...  STOP!!!", grid%xxbc(:,nc)
      
        print *,grid%varbc
        stop    
      endif 
    endif


  enddo

  call VecRestoreArrayF90(grid%icap, icap_p, ierr)
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
   
end subroutine pflow_Vadose_bcadj

end module Vadose_module
