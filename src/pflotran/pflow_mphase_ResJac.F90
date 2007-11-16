! introduced grid variables: e_total :: 1 dof
! translator_module                           c_total :: option%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (option%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

               
module MPHASE_module

  implicit none
  
  private

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscerror.h" 

! Cutoff parameters
!  real*8, parameter :: formeps   = 5.D-5
  real*8, parameter :: eps       = 1.D-5
  real*8, parameter :: floweps   = 1.D-24
!  real*8, parameter :: satcuteps = 1.D-5
  real*8, parameter :: zerocut =0.D0!1D-8
  real*8, parameter :: dfac = 1.D-8
!  real*8, parameter :: dfac = 1.D-2

  integer,save :: size_var_use 
  integer,save :: size_var_node
  real*8, allocatable,save :: Resold_AR(:,:), Resold_FL(:,:)
! Contributions to residual from accumlation/source/Reaction, flux(include diffusion)

  public MPHASEResidual, MPHASEJacobian, pflow_mphase_initaccum, &
         pflow_update_mphase,pflow_mphase_initadj, pflow_mphase_timecut,&
         pflow_mphase_setupini, MPhase_Update_Reason

  public :: createmphaseZeroArray
  
  integer, save :: n_zero_rows = 0
  integer, pointer, save :: zero_rows_local(:)  ! 1-based indexing
  integer, pointer, save :: zero_rows_local_ghosted(:) ! 0-based indexing
  
contains

subroutine pflow_mphase_timecut(realization)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module 
 
  implicit none

  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
 
  PetscScalar, pointer :: xx_p(:),yy_p(:)!,var_loc_p(:),iphase_loc_p(:)
  integer :: re,ierr
  integer :: dof_offset
  integer :: local_id
  !integer re0, ierr, index, iipha
  !real*8, pointer :: sat(:),xmol(:)

  grid => realization%grid
  option => realization%option
  field => realization%field  

  call VecGetArrayF90(field%xx, xx_p, ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr)
 !call VecGetArrayF90(field%var_loc, var_loc_p, ierr); 
 !call VecGetArrayF90(field%iphas, iphase_loc_p, ierr); 

  do local_id=1, grid%nlmax
    dof_offset=(local_id-1)*option%ndof
    do re = 1, option%ndof
!     if (option%snes_reason == -5) then 
!       xx_p(n0+re)= yy_p(n0+re)       !.5D0 * xx_p(n0+re) +.5D0 *yy_p(n0+re)
!     else
        xx_p(dof_offset+re)= yy_p(dof_offset+re)
!     endif
    enddo
  enddo 
  call VecRestoreArrayF90(field%xx, xx_p, ierr) 
  call VecRestoreArrayF90(field%yy, yy_p, ierr)
  
  !call VecCopy(field%xx,field%yy,ierr)
  !call pflow_mphase_initaccum(grid)
 
end subroutine pflow_mphase_timecut


subroutine pflow_mphase_setupini(realization)
 
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Region_module
  use Structured_Grid_module
  use Coupler_module
  use Condition_module
  
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: initial_condition
  
  PetscScalar, pointer :: xx_p(:), iphase_loc_p(:)
  integer local_id, ghosted_id, ibegin, iend, icell, ierr  
  
  grid => realization%grid
  option => realization%option
  field => realization%field
    
  size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
  size_var_node = (option%ndof + 1) * size_var_use
  
  allocate(Resold_AR(grid%nlmax,option%ndof))
  allocate(Resold_FL(grid%internal_connection_list%first%num_connections, &
                     option%ndof))
  allocate(option%delx(option%ndof,grid%ngmax))
  option%delx=0.D0
   
  call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  
  initial_condition => realization%initial_conditions%first
  
  do
  
    if (.not.associated(initial_condition)) exit
    
    do icell=1,initial_condition%region%num_cells
      local_id = initial_condition%region%cell_ids(icell)
      ghosted_id = grid%nL2G(local_id)
      iend = local_id*option%ndof
      ibegin = iend-option%ndof+1
      xx_p(ibegin:iend) = &
        initial_condition%condition%cur_value(1:option%ndof)
      iphase_loc_p(ghosted_id)=initial_condition%condition%iphase
    enddo
  
    initial_condition => initial_condition%next
  
  enddo
              
  call VecRestoreArrayF90(field%xx, xx_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p,ierr)

end subroutine pflow_mphase_setupini
  

subroutine MPhase_Update_Reason(reason,realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
    
  implicit none

  type(realization_type) :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
   
  integer, intent(out):: reason
  PetscScalar, pointer :: xx_p(:),var_loc_p(:),iphase_loc_p(:), yy_p(:) !,r_p(:)
  integer dof_offset, temp_reason
  integer :: ierr, iipha
  integer :: local_id, ghosted_id  

! real*8, pointer :: sat(:),xmol(:)
! real*8 rmax(option%ndof)

  grid => realization%grid
  option => realization%option
  field => realization%field
  
  reason=1
 ! call SNESComputeFunction(option%snes,field%xx,option%r,ierr)
 ! do n=1,option%ndof
 !  call VecStrideNorm(option%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if(reason>0)then
  call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p, ierr); 
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); 
  
  do local_id = 1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif  
     dof_offset=(local_id-1)* option%ndof
      !index=(n-1)*size_var_node
      !sat=>var_loc_p(index+2+1:index+2+option%nphase)
      !den=>var_loc_p(index+2+option%nphase+1:index+2+2*option%nphase)
    !xmol=>var_loc_p(index+2+7*option%nphase+1:index+2+7*option%nphase + option%nphase*option%nspec)    
      iipha=int(iphase_loc_p(local_id))
     !if(n==3583 .or. n==3587)
   !print *, 'update reson', grid%nlmax, n, iipha, xx_p(n0+1:n0+3)
   !if(xmol(4)>1.0) re=0; goto 1
   !if(xmol(4)<.0) re=0; goto 1
   !if(sat(2) < .0) re=0;goto 1
   ! if(sat(2) > 1.) re=0;goto 1

   if(dabs(xx_p(dof_offset + 1)- yy_p(dof_offset + 1))> (10.0D0 * option%dpmxe))then
     reason=0; print *,'huge change in p', xx_p(dof_offset + 1), yy_p(dof_offset + 1)
     exit
    endif


   if(dabs(xx_p(dof_offset + 2)- yy_p(dof_offset + 2))> (10.0D0 * option%dtmpmxe))then
     reason=0; print *,'huge change in T', xx_p(dof_offset + 2), yy_p(dof_offset + 2)
     exit
    endif
 
   select case(iipha)
   case (1)
     if(xx_p(dof_offset + 3) > 1.0D0)then
      reason=0; exit
!    goto 111
        endif
     if(xx_p(dof_offset + 3) < 0D0)then
      reason=0; exit
!    goto 111
       endif
     !if(xx_p(dof_offset + 3) > 1.0D0) xx_p(dof_offset + 3)=1.D0
     !if(xx_p(dof_offset + 3) < .0D0) xx_p(dof_offset + 3)=0.D0
    case (2)
     if(xx_p(dof_offset + 3) > 1.0D0)then
      reason=0; exit
!    goto 111
        endif
     if(xx_p(dof_offset + 3) < 0D-0)then
      reason=0; exit
!    goto 111
     endif
    case (3)
     if(xx_p(dof_offset + 3) > 1.D0)then
      reason=0; exit
!     goto 111
        endif
     if(xx_p(dof_offset + 3) < 0.)then
      reason=0; exit
!     goto 111
     endif
     !if(xx_p(dof_offset + 3) > 1.0D0) xx_p(dof_offset + 3)=1.D0
     !if(xx_p(dof_offset + 3) < .0D0) xx_p(dof_offset + 3)=0.D0
    end select  
  end do
  
!  do n = 1,grid%nlmax
!     n0=(n-1)* option%ndof
!      
!   if(dabs(xx_p(n0+1)-yy_p(n0+1))>1D6) then
!      re=0;exit
!   endif
   
!   if(dabs(xx_p(n0+2)-yy_p(n0+2))>1D1) then
!      re=0;exit
!   endif
  
  
!   enddo
   ! print *, 'update reason: ',option%myrank,grid%nlmax,n,re

  !   call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
   !print *,' update reason ba MPI', ierr
  if(reason<=0) print *,'Sat or Con out of Region at: ',local_id,iipha,xx_p(dof_offset+1:dof_offset+3)
    call VecRestoreArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%yy, yy_p, ierr)
    call VecRestoreArrayF90(field%var_loc, var_loc_p, ierr) 
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr) 
  endif
 ! print *,' update reason', option%myrank, re,n,grid%nlmax
  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  if(option%commsize >1)then
    temp_reason = reason
    call MPI_ALLREDUCE(temp_reason,reason,1, MPI_INTEGER,MPI_SUM, &
    PETSC_COMM_WORLD,ierr)
  !print *,' update reason re'
    !call MPI_BCAST(re0,1, MPI_INTEGER, 0,PETSC_COMM_WORLD,ierr)
  !print *,' update reason ca'
    if(reason<option%commsize) reason=0
  endif
  
  if(reason<=0 .and.option%myrank ==0) print *,'Sat or Con out of Region', reason
  
end subroutine MPhase_Update_Reason




subroutine MPHASERes_ARCont(node_no, var_node,por,vol,rock_dencpr, option, &
                              Res_AR,ireac,ierr)
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  integer :: node_no
  integer, optional:: ireac,ierr
  real*8, target :: var_node(1:size_var_use)
  real*8 :: Res_AR(1:option%ndof) 
  real*8 :: vol,por,rock_dencpr
     
  real*8, pointer :: temp, pre_ref   ! 1 dof
  real*8, pointer :: sat(:), density(:), amw(:), h(:), u(:), pc(:), kvr(:) ! nphase dof
  real*8, pointer :: xmol(:), diff(:)            ! nphase*nspec
  
  integer :: ibase, m,np, iireac=1
  real*8 :: pvol,mol(option%nspec),eng
  
  if(present(ireac)) iireac=ireac
  pvol=vol*por
  
  ibase=1;                 temp=>var_node(ibase)
  ibase=ibase+1;           pre_ref=>var_node(ibase)
  ibase=ibase+1;           sat=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; density=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; amw=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; h=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; u=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; pc=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; kvr=>var_node(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; xmol=>var_node(ibase:ibase+option%nphase*option%nspec-1)
  ibase=ibase+option%nphase*option%nspec; 
  diff=>var_node(ibase:ibase+option%nphase*option%nspec-1)

  !sumation of component
  mol=0.D0; eng=0.D0
  do np = 1, option%nphase
    do m=1, option%nspec  
      mol(m) = mol(m) + sat(np)*density(np)*xmol(m + (np-1)*option%nspec)
    enddo
    eng = eng + density(np)*u(np) *sat(np)
  enddo

  mol = mol * pvol
  eng = eng * pvol + (1.D0 - por)* vol * rock_dencpr * temp 
  
! Reaction terms here
  if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!   H2O
    mol(1)= mol(1) - option%dt * option%rtot(node_no,1)
!   CO2
    mol(2)= mol(2) - option%dt * option%rtot(node_no,2)
!   should include related energy change here
  endif
  Res_AR(1:option%ndof-1)=mol(:)
  Res_AR(option%ndof)=eng

  nullify(temp, pre_ref, sat, density, amw, h,u, pc,kvr,xmol,diff)
  
end subroutine  MPHASERes_ARCont


subroutine MPHASERes_FLCont(nconn_no,area, &
                             var_node1,por1,tor1,sir1,dd1,perm1,Dk1,&
                             var_node2,por2,tor2,sir2,dd2,perm2,Dk2,dist_gravity,upweight,&
                             option, vv_darcy,Res_FL)
  
  use Option_module                             

  implicit none
  
  integer nconn_no
  type(option_type) :: option
  real*8 :: sir1(1:option%nphase),sir2(1:option%nphase)
  real*8, target :: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8, target :: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8 :: por1,por2,tor1,tor2,perm1,perm2,Dk1,Dk2,dd1,dd2
  real*8 :: vv_darcy(option%nphase),area
  real*8 :: Res_FL(1:option%ndof) 
  
  real*8 :: dist_gravity ! distance along gravity vector
     
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  real*8, pointer :: xmol1(:), diff1(:)            ! 
  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  real*8, pointer :: xmol2(:), diff2(:)    
  
  integer :: ibase, m,np, ind
  real*8 :: fluxm(option%nspec),fluxe, v_darcy,q
  real*8 :: uh,uxmol(1:option%nspec), ukvr,difff,diffdp, DK,Dq
  real*8 :: upweight,density_ave,cond, gravity, dphi
  
!  m1=option%nd1(nc); n1 = grid%nG2L(m1) ! = zero for ghost nodes 
!  print *,'in FLcont'
  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
               
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
               
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+option%nphase-1)
               sat2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; density1=>var_node1(ibase:ibase+option%nphase-1)
                           density2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; amw1=>var_node1(ibase:ibase+option%nphase-1)
                           amw2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; h1=>var_node1(ibase:ibase+option%nphase-1)
                           h2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; u1=>var_node1(ibase:ibase+option%nphase-1)
                           u2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; pc1=>var_node1(ibase:ibase+option%nphase-1)
                           pc2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; kvr1=>var_node1(ibase:ibase+option%nphase-1)
                           kvr2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; xmol1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)
                           xmol2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)
  ibase=ibase+option%nphase*option%nspec;
               diff1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)    
               diff2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)

  !print *,' FLcont got pointers' ,var_node1,var_node2,sir1,sir2
  !print *,' tmp=',temp1,temp2
  !print *,'diff=',diff1,diff2
   
  Dq=(perm1 * perm2)/(dd1*perm2 + dd2*perm1)
  diffdp= (por1 *tor1 * por2*tor2) / (dd2*por1*tor1 + dd1*por2*tor2)*area
  
  fluxm=0.D0
  fluxe=0.D0
  vv_darcy=0.D0  
  
  do np=1, option%nphase

!   Flow term
    if ((sat1(np) > sir1(np)) .or. (sat2(np) > sir2(np)))then
    
      upweight= dd2/(dd1+dd2)
      if(sat1(np) <eps) then 
        upweight=0.d0
      else if(sat2(np) <eps) then 
        upweight=1.d0
      endif
      density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
      
      gravity = (upweight*density1(np)*amw1(np) + &
      (1.D0-upweight)*density2(np)*amw2(np)) &
      * option%gravity * dist_gravity
          
      dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
!     print *,'FLcont  dp',dphi
  !   note uxmol only contains one phase xmol
  !   removed iupstream for preserving derivative wrt pressure
!     if(option%iupstream(0,1) ==0)then
        if(dphi>=0.D0)then
          ukvr=kvr1(np)
          uh=h1(np)
          uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
        else
          ukvr=kvr2(np)
          uh=h2(np)
          uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
        endif      
!     else
!       if(option%iupstream(nconn_no, np) >=0 )then
!         ukvr=kvr1(np)
!         uh=h1(np)
!         uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
!       else
!         ukvr=kvr2(np)
!         uh=h2(np)
!         uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
!       endif      
!     endif  
     
     !print *,'FLcont  uxmol',uxmol
      if(ukvr>floweps)then
        v_darcy= Dq * ukvr * dphi
        !option%vvl_loc(nconn_no) = v_darcy
        vv_darcy(np)=v_darcy
        q=v_darcy * area
        do m=1, option%nspec 
          fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
        enddo
        fluxe = fluxe + q*density_ave*uh 
      endif
    endif 
 !  print *,' FLcont end flow',np
!   Diffusion term   
!   Note : average rule may not be correct  
    if ((sat1(np) > eps) .and. (sat2(np) > eps))then
      difff =diffdp * 0.25D0*(sat1(np)+sat2(np))*(density1(np)+density2(np))
      do m = 1, option%nspec
        ind=m+(np-1)*option%nspec
        fluxm(m) = fluxm(m) + difff * .5D0 * &
        (diff1(ind) + diff2(ind))*(xmol1(ind)-xmol2(ind))
      enddo  
    endif 
  enddo
   
! conduction term
        
  Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)
  cond=Dk*area*(temp1-temp2) 
  fluxe=fluxe + cond
   !print *,' FLcont heat cond', Dk, cond
  Res_FL(1:option%ndof-1)=fluxm(:) * option%dt
  Res_FL(option%ndof)=fluxe * option%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 ! 2 R = R - Res_FL
 !print *,'end FLcont'

  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)

end subroutine MPHASERes_FLCont


subroutine MPHASERes_FLBCCont(nbc_no,ibndtype,area, &
            var_node1,var_node2,por2,tor2,sir2,dd1,perm2,Dk2,dist_gravity,&
            option,field,vv_darcy,Res_FL)

  use Option_module
  use Field_module
  
  implicit none
  
  integer nbc_no, ibndtype
  type(option_type) :: option
  type(field_type) :: field
  real*8 :: dd1, sir2(1:option%nphase)
  real*8, target :: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8, target :: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8 :: por2,perm2,Dk2,tor2
  real*8 :: vv_darcy(option%nphase), area
  real*8 :: Res_FL(1:option%ndof) 
  
  real*8 :: dist_gravity   ! distance along gravity vector
     
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  real*8, pointer :: xmol1(:), diff1(:)            ! 
  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  real*8, pointer :: xmol2(:), diff2(:)    
  
  integer :: ibase, m,np, ind, ibc,j
  real*8 :: fluxm(option%nspec),fluxe, v_darcy,q
  real*8 :: uh,uxmol(1:option%nspec), ukvr,diff,diffdp, DK,Dq
  real*8 :: upweight,density_ave,cond,gravity, dphi

  
  ibase=1;                 temp1=>var_node1(ibase)
                           temp2=>var_node2(ibase)
  ibase=ibase+1;           pre_ref1=>var_node1(ibase)
                           pre_ref2=>var_node2(ibase)
  ibase=ibase+1;           sat1=>var_node1(ibase:ibase+option%nphase-1)
               sat2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; density1=>var_node1(ibase:ibase+option%nphase-1)
                           density2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; amw1=>var_node1(ibase:ibase+option%nphase-1)
                           amw2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; h1=>var_node1(ibase:ibase+option%nphase-1)
                           h2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; u1=>var_node1(ibase:ibase+option%nphase-1)
                           u2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; pc1=>var_node1(ibase:ibase+option%nphase-1)
                           pc2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; kvr1=>var_node1(ibase:ibase+option%nphase-1)
                           kvr2=>var_node2(ibase:ibase+option%nphase-1)
  ibase=ibase+option%nphase; xmol1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)
                           xmol2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)
  ibase=ibase+option%nphase*option%nspec;
               diff1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)    
               diff2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)


  fluxm=0.D0
  fluxe=0.D0
  vv_darcy=0.D0 
   
  select case (ibndtype)
  
    case(1)
    
      Dq= perm2 / dd1
      diffdp = por2*tor2/dd1*area
    ! Flow term
      do np =1,option%nphase
        if ((sat1(np) > sir2(np)) .or. (sat2(np) > sir2(np)))then
    
          upweight=1.D0
          if(sat1(np) <eps) then 
            upweight=0.d0
          else if(sat2(np) <eps) then 
            upweight=1.d0
          endif
          density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
      
          gravity = (upweight*density1(np)*amw1(np) + &
          (1.D0-upweight)*density2(np)*amw2(np)) &
          * option%gravity * dist_gravity
          
          dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
    
!         if(option%iupstream(0,1) ==0)then
            if(dphi>=0.D0)then
              ukvr=kvr1(np)
              uh=h1(np)
              uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
           else
              ukvr=kvr2(np)
              uh=h2(np)
              uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
            endif      
!         else
!          if(option%iupstreambc(nbc_no, np) >=0 )then
!             ukvr=kvr1(np)
!             uh=h1(np)
!             uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
!           else
!             ukvr=kvr2(np)
!             uh=h2(np)
!             uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
!           endif      
!         endif  


          if(ukvr>floweps)then
            v_darcy= Dq * ukvr * dphi
            !option%vvl_loc(nbc_no) = v_darcy
            vv_darcy(np)=v_darcy
     
            q=v_darcy * area
          
            do m=1, option%nspec 
              fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
            enddo
            fluxe = fluxe + q*density_ave*uh 
          endif
        endif
        
!       Diffusion term   
!       Note : average rule may not be correct  
        if ((sat1(np) > eps) .and. (sat2(np) > eps))then
          diff =diffdp * 0.25D0*(sat1(np)+sat2(np))*(density1(np)+density2(np))
          do m = 1, option%nspec
            ind=m+(np-1)*option%nspec
            fluxm(m) = fluxm(m) + diff * diff2(ind)&
            *( xmol1(ind)-xmol2(ind))
          enddo  
        endif
      enddo
      
!     conduction term

      Dk =  Dk2 / dd1
      cond=Dk*area*(temp1-temp2) 
      fluxe=fluxe + cond
 
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%ndof)=fluxe * option%dt

    case(2)
    
      if((dabs(field%velocitybc(1,nbc_no))+ &
      dabs(field%velocitybc(2,nbc_no)))>floweps)then
!       print *, 'FlowBC :',nbc_no,field%velocitybc(1,nbc_no),field%velocitybc(2,nbc_no)
        do j=1,option%nphase
 !        fluxm=0.D0; fluxe=0.D0
          v_darcy = field%velocitybc(j,nbc_no)
          vv_darcy(j) = field%velocitybc(j,nbc_no)
!         option%vvbc(j+(nc-2)*option%nphase)= field%velocitybc(j,nc)
      !   note different from 2 phase version

          if(v_darcy >0.d0)then 
            q = v_darcy * density1(j) * area
            !q = 0.d0
            !flux = flux - q
            fluxe = fluxe - q  * h1(j) 
            do m=1, option%nspec
              fluxm(m)=fluxm(m) + q * xmol1(m + (j-1)*option%nspec)
            enddo 
          else 
            q =  v_darcy * density2(j) * area   
            fluxe = fluxe - q  * h2(j) 
            do m=1, option%nspec
              fluxm(m)=fluxm(m) + q * xmol2(m + (j-1)*option%nspec)
            enddo 
          endif 
        enddo
      endif
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%ndof)=fluxe * option%dt

    case(3)
   
      Dq= perm2 / dd1
      diffdp = por2*tor2/dd1*area
      
    ! Flow term
      do np =1,option%nphase
        if ((sat1(np) > sir2(np)) .or. (sat2(np) > sir2(np)))then
          upweight=1.D0
          if(sat1(np) <eps) then 
            upweight=0.d0
          else if(sat2(np) <eps) then 
            upweight=1.d0
          endif
          density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  
    
          gravity = (upweight*density1(np)*amw1(np) + &
          (1.D0-upweight)*density2(np)*amw2(np)) &
              * option%gravity * dist_gravity
        
          dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
    
!         if(option%iupstream(0,1) ==0)then
             if(dphi>=0.D0)then
               ukvr=kvr1(np)
               uh=h1(np)
               uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
            else
               ukvr=kvr2(np)
               uh=h2(np)
               uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
             endif      
!          else
!           if(option%iupstreambc(nbc_no, np) >=0 )then
!              ukvr=kvr1(np)
!              uh=h1(np)
!              uxmol(1:option%nspec)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
!            else
!              ukvr=kvr2(np)
!              uh=h2(np)
!              uxmol(1:option%nspec)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
!            endif      
!          endif  
 
         
          if(ukvr>floweps)then
            v_darcy= Dq * ukvr * dphi
               !option%vvl_loc(nbc_no) = v_darcy
            vv_darcy(np)=v_darcy
           
            q=v_darcy * area
                
            do m=1, option%nspec 
              fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
            enddo
            fluxe = fluxe + q*density_ave*uh 
          endif
        endif 
      enddo
   
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%ndof)=fluxe * option%dt
    
    case(4)
          
      Dk =  Dk2 / dd1
      cond = Dk*area*(temp1-temp2) 
      fluxe=fluxe + cond
     
      Res_FL(1:option%nspec)= 0.D0
      Res_FL(option%ndof)=fluxe * option%dt

  end select
  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)
  
end  subroutine MPHASERes_FLBCCont 


subroutine MPHASEResidual(snes,xx,r,realization,ierr)

  use water_eos_module
  use co2eos_module
  use translator_mph_module
  use span_wagner_module

  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module 
  use Field_module

  implicit none
  
#include "definitions.h"  

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  integer :: ierr, ierr0
  integer :: nr
  integer :: i, jn
  integer :: ip2, p1, p2
  integer :: ithrm1, ithrm2
  integer :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  
  PetscScalar, pointer :: accum_p(:)

  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), var_loc_p(:) 
                          

  PetscScalar, pointer :: iphase_loc_p(:),icap_loc_p(:), ithrm_loc_p(:)

  integer :: iicap,iiphase, index_var_begin, index_var_end,iicap1,iicap2,np

  real*8 :: dd1, dd2, eng, &
!           eengl,eengg, &
!           fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl,&
!           fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q,  &
!           v_darcy,hflx,
            pvoldt, voldt, accum, pvol
  real*8 :: dd, f1, f2, ff, perm1, perm2
! real*8 :: Dphi,D0, por1, por2,density_ave
! real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
! real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol,dif(realization%option%nphase)
  real*8 :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1!, qqsrc
! real*8 :: cw1,cw2, xxlw,xxla,xxgw,xxga
  real*8 :: cw
! real*8 :: upweight
! real*8 :: ukvr,uhh,uconc
  real*8 :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp, rho, visc
  real*8 :: Res(realization%option%ndof), vv_darcy(realization%option%nphase)
 
! real*8 :: cond, den,
  PetscViewer :: viewer

  ! Connection variables
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  logical :: enthalpy_flag
  integer :: iconn
  integer :: sum_connection
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity, upweight

  grid => realization%grid
  option => realization%option
  field => realization%field
 
#if 0 
 ! only vvl_loc is used, and it is commented out   
  field%vvlbc=0.D0
  field%vvgbc=0.D0
  field%vvl_loc=0.D0
  field%vvg_loc=0.D0
#endif  

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  
  ierr = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
      
!   insure zero liquid sat not passed to ptran (no effect on pflow)
    if(xx_p((local_id-1)*option%ndof+3) < 0.D0)xx_p((local_id-1)*option%ndof+3) = zerocut
    if(xx_p((local_id-1)*option%ndof+3) > 1.D0)xx_p((local_id-1)*option%ndof+3) = 1.D0 - zerocut
    
!   check if p,T within range of table  
    if(xx_p((local_id-1)*option%ndof+1)< p0_tab*1D6 &
      .or. xx_p((local_id-1)*option%ndof+1)>(ntab_p*dp_tab + p0_tab)*1D6)then
       ierr=-1; exit  
     endif  
    if(xx_p((local_id-1)*option%ndof+2)< t0_tab -273.15D0 &
      .or. xx_p((local_id-1)*option%ndof+2)>ntab_t*dt_tab + t0_tab-273.15D0)then
      ierr=-1; exit
    endif
  enddo
  
  ierr0 = 0
  if(option%commsize >1)then
    call MPI_ALLREDUCE(ierr, ierr0,1, MPI_INTEGER,MPI_SUM, PETSC_COMM_WORLD,ierr)
    if(ierr0 < 0) then
      ierr=-1      
    endif
  endif
  
  if(ierr<0)then
    ierr = PETSC_ERR_ARG_DOMAIN
    if (option%myrank==0) print *,'table out of range: ',ierr0
    return
  endif 


! allow phase change for first 3 newton iterations except zeroth iteration
  if(option%iphch>0 .and. option%iphch<=3)then
!  if(option%iphch<=3)then
    call Translator_MPhase_Switching(xx,realization,0,ierr)   
  endif  
  option%iphch=option%iphch+1
   
  call VecRestoreArrayF90(xx, xx_p, ierr)

  call GridGlobalToLocal(grid,xx,field%xx_loc,NDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)                          

  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); CHKERRQ(ierr)

! there is potential possiblity that the pertubation of p may change the direction of pflow.
! once that happens, code may crash, namely wrong derive. 
  do ghosted_id = 1, grid%ngmax 

    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
      
    iiphase=int(iphase_loc_p(ghosted_id))

!#if 0
    option%delx(1,ghosted_id) = xx_loc_p((ghosted_id-1)*option%ndof+1)*dfac * 1.D-3
    option%delx(2,ghosted_id) = xx_loc_p((ghosted_id-1)*option%ndof+2)*dfac

    select case (iiphase)
      case (1)
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) < .8)then
          option%delx(3,ghosted_id) =  dfac*xx_loc_p((ghosted_id-1)*option%ndof+3)
        else
          option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%ndof+3) 
        endif
        if(option%delx(3,ghosted_id) < 1D-9 .and. option%delx(3,ghosted_id)>=0.D0)option%delx(3,ghosted_id) =1D-9
        if(option%delx(3,ghosted_id) >-1D-9 .and. option%delx(3,ghosted_id)<0.D0)option%delx(3,ghosted_id) =-1D-9
      case(2)  
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) <0.8)then
          option%delx(3,ghosted_id) =  dfac*xx_loc_p((ghosted_id-1)*option%ndof+3) 
        else
          option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%ndof+3) 
        endif 
        if(option%delx(3,ghosted_id) < 1D-9 .and. option%delx(3,ghosted_id)>=0.D0)option%delx(3,ghosted_id) =1D-9
        if(option%delx(3,ghosted_id) >-1D-9 .and. option%delx(3,ghosted_id)<0.D0)option%delx(3,ghosted_id) =-1D-9
      case(3)
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) <=0.9)then
          option%delx(3,ghosted_id) = dfac*xx_loc_p((ghosted_id-1)*option%ndof+3) 
        else
          option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%ndof+3) 
        endif 
        
        if(option%delx(3,ghosted_id) < 1D-9 .and. option%delx(3,ghosted_id)>=0.D0)option%delx(3,ghosted_id) = 1D-9
        if(option%delx(3,ghosted_id) >-1D-9 .and. option%delx(3,ghosted_id)<0.D0)option%delx(3,ghosted_id) =-1D-9
        
        if((option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%ndof+3))>1.D0)then
          option%delx(3,ghosted_id) = (1.D0-xx_loc_p((ghosted_id-1)*option%ndof+3))/1D5
        endif
        if((option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%ndof+3))<0.D0)then
          option%delx(3,ghosted_id) = xx_loc_p((ghosted_id-1)*option%ndof+3)/1D5
        endif
    end select
!#endif

#if 0
    option%delx(1,ng) = 1D-1
    option%delx(2,ng) = 1D-8

    select case (iiphase)
      case (1)
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) < 5D-5)then
          option%delx(3,ghosted_id) = 1D-8
        else
          option%delx(3,ghosted_id) = -1D-8 
        endif
      case(2)  
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) <0.9995)then
          option%delx(3,ghosted_id) =  1D8
        else
          option%delx(3,ghosted_id) = -1D-8 
        endif 
      case(3)
        if(xx_loc_p((ghosted_id-1)*option%ndof+3) <=0.9)then
          option%delx(3,ghosted_id) = 1D-10 
        else
          option%delx(3,ghosted_id) = -1D-10 
        endif 
        
        
        if((option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%ndof+3))>1.D0)then
          option%delx(3,ghosted_id) = (1.D0-xx_loc_p((ghosted_id-1)*option%ndof+3))*1D-6
        endif
        if((option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%ndof+3))<0.D0)then
          option%delx(3,ghosted_id) = xx_loc_p((ghosted_id-1)*option%ndof+3)*1D-6
        endif
    end select
#endif

#if 0
    option%delx(1,ghosted_id) = 1.d-1 ! pressure increment
    option%delx(2,ghosted_id) = 1.d-7 ! temperature increment

    if(iiphase==2)then
       if(xx_loc_p((ng-1)*option%ndof+3) <=0.9995)then
        option%delx(3,ghosted_id) =  1.d-9
      else
        option%delx(3,ghosted_id) = -1.d-8
      endif
    endif

    if(iiphase==1)then
       if(xx_loc_p((ng-1)*option%ndof+3) <=5D-4)then
        option%delx(3,ghosted_id) =  1.d-8
      else
        option%delx(3,ghosted_id) = -1.d-8
      endif
    endif
       
    if (iiphase == 3) then
      option%delx(3,ghosted_id) = -1.d-8
      if (xx_loc_p((ghosted_id-1)*option%ndof+3) <= 0.01) option%delx(3,ghosted_id) = 1.d-8
    endif
#endif
  enddo
  
  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%var_loc,var_loc_p,ierr)
  
  
 ! call VecGetArrayF90(field%ithrm,ithrm_loc_p,ierr)
!------------------------------------------------------ 





!-----  phase properities ---- last time step---
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif  
    
    jn = 1 + (local_id-1)*option%nphase
    iicap = int(icap_loc_p(ghosted_id))
    iiphase = iphase_loc_p(ghosted_id)
    !*****************
    dif(1)= option%difaq
    dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
 ! print *, n,iicap,xx_p((n-1)*option%ndof+1:n*option%ndof) 
  !*******************************************
    call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase,&
        option%scale,option%nphase,option%nspec,iicap, dif,&
    var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),&
    option%itable,option%m_nacl,ierr,field%xxphi_co2(local_id), field%dden_co2(local_id))



    if (option%ideriv .eq. 1) then
      call pri_var_trans_mph_winc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),&
        option%delx(1:option%ndof,ghosted_id), iiphase,&
        option%scale,option%nphase,option%nspec, iicap, dif,&
        var_loc_p((ghosted_id-1)*size_var_node+size_var_use+1:ghosted_id*size_var_node),&
      option%itable,option%m_nacl,ierr)
    endif
  
!  print *,'var_loc_p',n,iicap,iiphase, var_loc_p((n-1)*size_var_node+1:n*size_var_node)              
!   if(n < 5) print *,'pflow_2ph: ',n,option%ideriv,field%xxphi_co2(n)
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecRestoreArrayF90(field%var_loc,var_loc_p,ierr)
  ! call VecRestoreArrayF90(field%iphas_loce,iphase_loc_p,ierr)
  

  call GridLocalToLocal(grid,field%var_loc,field%var_loc,VARDOF)
  call GridLocalToLocal(grid,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)


! End distribute data 
! now assign access pointer to local variables
  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%accum, accum_p, ierr)
! call VecGetArrayF90(field%yy, yy_p, ierr)
 

  ! notice:: here we assume porosity is constant
 
  call VecGetArrayF90(field%var_loc,var_loc_p,ierr)
  call VecGetArrayF90(field%yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%vl, vl_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'


  if (option%rk > 0.d0) then
    call VecGetArrayF90(field%phis,phis_p,ierr)
  endif

  Resold_AR=0.D0; ResOld_FL=0.D0

!--------------------------------------------------------------------------
! Calculate accumulation term for interior and exterior nodes.
!--------------------------------------------------------------------------
! print *,option%rtot
  
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...

    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    
    p1 = 1 + (local_id-1)*option%ndof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    
    pvol = volume_p(local_id)*porosity_loc_p(ghosted_id)
    voldt = volume_p(local_id) / option%dt
    pvoldt = porosity_loc_p(ghosted_id) * voldt
    iiphase = iphase_loc_p(ghosted_id)
    i = ithrm_loc_p(ghosted_id)

    accum = 0.d0
    call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin: index_var_end),&
    porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(i), option, Res, 1,ierr)
   
    r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
    Resold_AR(local_id,1:option%ndof)= Res(1:option%ndof) 
  end do

!************************************************************************
! add source/sink terms
#if 0 
  do nr = 1, option%nblksrc
      
    kk1 = option%k1src(nr) - option%nzs
    kk2 = option%k2src(nr) - option%nzs
    jj1 = option%j1src(nr) - option%nys
    jj2 = option%j2src(nr) - option%nys
    ii1 = option%i1src(nr) - option%nxs
    ii2 = option%i2src(nr) - option%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(option%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(option%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(option%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, option%ntimsrc
      if (option%timesrc(i,nr) == option%t) then
        tsrc1 = option%tempsrc(i,nr)
        qsrc1 = option%qsrc(i,nr)
        csrc1 = option%csrc(i,nr)
        hsrc1 = option%hsrc(i,nr)
        goto 10
      else if (option%timesrc(i,nr) > option%t) then
        ff = option%timesrc(i,nr)-option%timesrc(i-1,nr)
        f1 = (option%t - option%timesrc(i-1,nr))/ff
        f2 = (option%timesrc(i,nr)-option%t)/ff
        tsrc1 = f1*option%tempsrc(i,nr) + f2*option%tempsrc(i-1,nr)
        qsrc1 = f1*option%qsrc(i,nr) + f2*option%qsrc(i-1,nr)
        csrc1 = f1*option%csrc(i,nr) + f2*option%csrc(i-1,nr)
        hsrc1 = f1*option%hsrc(i,nr) + f2*option%hsrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
   !print *,'pflow2ph : ', option%myrank,i,option%timesrc(i,nr), &
   !option%timesrc(i-1,nr),option%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
  
  if(dabs(hsrc1)>1D-20)then 
    do kk = kk1, kk2
      do jj = jj1, jj2
        do ii = ii1, ii2
          n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
          r_p(n*option%ndof) = r_p(n*option%ndof) - hsrc1 * option%dt   
        enddo
      enddo
    enddo
  endif         
     
    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
            ng = grid%nL2G(n)

            call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+2),&
      dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
            
            r_p((n-1)*option%ndof + option%jh2o) = r_p((n-1)*option%ndof +option%jh2o) - qsrc1 *option%dt
            r_p(n*option%ndof) = r_p(n*option%ndof) - qsrc1*enth_src_h2o*option%dt
            Resold_AR(n,option%jh2o)= Resold_AR(n,option%jh2o) - qsrc1*option%dt
      Resold_AR(n,option%ndof)= Resold_AR(n,option%ndof) - qsrc1 * enth_src_h2o*option%dt
      
      
      !           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*option%fmwh2o, &
!           qsrc1
          enddo
        enddo
      enddo
    endif  
    
    if (csrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*option%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            rho = var_loc_p((ghosted_id-1)*size_var_node+4+option%nphase)*option%fmwco2 
            call co2_span_wagner(var_loc_p((ghosted_id-1)*size_var_node+2)*1.D-6,&
                  tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                        eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)

         !  units: rho [kg/m^3]; csrc1 [kmol/s]

            enth_src_co2 = enth_src_co2 * option%fmwco2
           
            r_p((n-1)*option%ndof + option%jco2) = r_p((n-1)*option%ndof + option%jco2) - csrc1*option%dt
            r_p(n*option%ndof) = r_p(n*option%ndof) - csrc1 * enth_src_co2 *option%dt
            Resold_AR(n,option%jco2)= Resold_AR(n,option%jco2) - csrc1*option%dt
      Resold_AR(n,option%ndof)= Resold_AR(n,option%ndof) - csrc1 * enth_src_co2*option%dt
       !r_p(s1) = r_p(s1) - csrc1

        !   print *,'pflow2ph_co2: ',option%myrank,nr,n,ng,tsrc1,rho,option%fmwco2,csrc1
          enddo
        enddo
      enddo
    endif
  
  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
      
    !  do kk = kk1, kk2
     !   do jj = jj1, jj2
       !   do ii = ii1, ii2
          !    n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
           !   ng = grid%nL2G(n)
           !   p1 = 1+(n-1)*option%ndof
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
 ! print *,'finished source/sink term'
#endif  

!*********************************************************************


 
! stop
!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------
! option%iupstream = 0; option%iupstreambc = 0

  ! loop over internal connections
  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      ghosted_id_up = cur_connection_object%id_up(iconn)
      ghosted_id_dn = cur_connection_object%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or.  &
            field%imat(ghosted_id_dn) <= 0) cycle
      endif

      p1 = 1 + (local_id_up-1)*option%ndof 
      p2 = 1 + (local_id_dn-1)*option%ndof

      fraction_upwind = cur_connection_object%dist(-1,iconn)
      distance = cur_connection_object%dist(0,iconn)
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = cur_connection_object%dist(3,iconn)*distance
      dd1 = distance*fraction_upwind
      dd2 = distance-dd1 ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd2/(dd1+dd2)
        
      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(3,iconn))

      perm2 = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(3,iconn))

      ithrm1 = ithrm_loc_p(ghosted_id_up)
      ithrm2 = ithrm_loc_p(ghosted_id_dn)
      iicap1=int(icap_loc_p(ghosted_id_up))
      iicap2=int(icap_loc_p(ghosted_id_dn))
   
      D1 = option%ckwet(ithrm1)
      D2 = option%ckwet(ithrm2)

      call MPHASERes_FLCont(iconn,cur_connection_object%area(iconn), &
                            var_loc_p((ghosted_id_up-1)*size_var_node+1: &
                                      (ghosted_id_up-1)*size_var_node+size_var_use),&
                            porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((ghosted_id_dn-1)*size_var_node+1: &
                                      (ghosted_id_dn-1)*size_var_node+size_var_use),&
                            porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2,&
                            distance_gravity,upweight,option,vv_darcy,Res)
      cur_connection_object%velocity(1,iconn) = vv_darcy(1) ! liquid
      cur_connection_object%velocity(2,iconn) = vv_darcy(2) ! gas
    
!   do np = 1, option%nphase
!     if(vv_darcy(np)>=0.D0) option%iupstream(nc,np) = 1
!     if(vv_darcy(np)<0.D0) option%iupstream(nc,np) = -1
!   enddo

      if (local_id_up > 0) then               ! If the upstream node is not a ghost node...
        do np =1, option%nphase 
          vl_p(np+(0)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                              vv_darcy(np)*cur_connection_object%dist(1,iconn) 
          vl_p(np+(1)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                              vv_darcy(np)*cur_connection_object%dist(2,iconn) 
          vl_p(np+(2)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                              vv_darcy(np)*cur_connection_object%dist(3,iconn) 
        enddo
      endif
     
      Resold_FL(iconn,1:option%ndof) = Res(1:option%ndof) 
    
      if(local_id_up>0)then
        r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
      endif
   
      if(local_id_dn>0)then
        r_p(p2:p2+option%ndof-1) = r_p(p2:p2+option%ndof-1) - Res(1:option%ndof)
      endif

    enddo
  
    cur_connection_object => cur_connection_object%next

  end do

 
!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', option%myrank, option%sgbc    

 
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection
    
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if(ghosted_id<=0)then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      end if

      p1 = 1 + (local_id-1) * option%ndof

      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_object%dist(3,iconn))
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = cur_connection_object%dist(3,iconn) * &
                         cur_connection_object%dist(0,iconn)

      select case(boundary_condition%condition%itype(1))
          
        case(2)
        ! solve for pb from Darcy's law given qb /= 0
          field%xxbc(:,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+1: ghosted_id*option%ndof)
          field%iphasebc(sum_connection) = int(iphase_loc_p(ghosted_id))
        case(3) 
       !  field%xxbc((sum_connection-1)*option%ndof+1)=option%pressurebc(2,ibc)
          field%xxbc(2:option%ndof,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+2: ghosted_id*option%ndof)
          field%iphasebc(sum_connection)=int(iphase_loc_p(ghosted_id))
        case(4)
           field%xxbc(1,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+1)
           field%xxbc(3:option%ndof,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+3: ghosted_id*option%ndof)    
           field%iphasebc(sum_connection)=int(iphase_loc_p(ghosted_id))

      end select
    
      iicap=int(icap_loc_p(ghosted_id))  
       
    !*****************
      dif(1)= option%difaq
      dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    !*******************************************

  
      call pri_var_trans_mph_ninc(field%xxbc(:,sum_connection),field%iphasebc(sum_connection), &
                                  option%scale,option%nphase,option%nspec,iicap, &
                                  dif,field%varbc(1:size_var_use),option%itable, &
                                  option%m_nacl,ierr,field%xxphi_co2_bc(sum_connection),cw)
     
      call MPHASERes_FLBCCont(sum_connection,boundary_condition%condition%itype(1), &
                              cur_connection_object%area(iconn), &
                              field%varbc(1:size_var_use), &
                              var_loc_p((ghosted_id-1)*size_var_node+1: &
                                        (ghosted_id-1)*size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id),tor_loc_p(ghosted_id), &
                              option%sir(1:option%nphase,iicap), &
                              cur_connection_object%dist(0,iconn),perm1,D2, &
                              distance_gravity,option,field,vv_darcy,Res)
                              
      cur_connection_object%velocity(1,iconn) = vv_darcy(1)  ! liquid
      cur_connection_object%velocity(2,iconn) = vv_darcy(2)  ! gas

      r_p(p1:p1-1+option%ndof)= r_p(p1:p1-1+option%ndof) - Res(1:option%ndof)
      ResOld_AR(local_id,1:option%ndof) = ResOld_AR(local_id,1:option%ndof) - Res(1:option%ndof)
    enddo
    boundary_condition => boundary_condition%next
  enddo
 

!adjust residual to R/dt

  select case (option%idt_switch) 
    case(1) 
      r_p(:) = r_p(:)/option%dt
    case(-1)
      if(option%dt>1.D0) r_p(:) = r_p(:)/option%dt
  end select    
  
  do local_id = 1, grid%nlmax
    p1 = 1 + (local_id-1)*option%ndof
    if(volume_p(local_id)>1.D0) r_p (p1:p1+2)=r_p(p1:p1+2)/volume_p(local_id)
  enddo  

  ! for inactive regions
  if (option%use_isoth==PETSC_TRUE) then
    print *, 'option%use_isoth needs to be verified in pflow_mphase_ResJac.F90'
    stop
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif
      p1 = 3 + (local_id-1)*option%ndof
      r_p(p1)=xx_loc_p(2 + (ghosted_id-1)*option%ndof)-yy_p(p1-1)
    enddo
  endif

  if (n_zero_rows > 0) then
    do i=1,n_zero_rows
      r_p(zero_rows_local(i)) = 0.
    enddo
  endif
   
!print *,'res closing pointer'
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr)
  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%accum, accum_p, ierr)
  call VecRestoreArrayF90(field%var_loc,var_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%vl, vl_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif

!#define DEBUG_GEH
!#define DEBUG_GEH_ALL
#ifdef DEBUG_GEH 
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'residual.out',viewer,ierr)
 call VecView(r,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'xx.out',viewer,ierr)
 call VecView(xx,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine MPHASEResidual
                
! --------------------------------------------------------------------- 

subroutine MPHASEJacobian(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module
  use co2eos_module
  use translator_mph_module, only : pri_var_trans_mph_ninc, pri_var_trans_mph_winc
  use span_wagner_module

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Coupler_module
  
  implicit none

#include "definitions.h"  

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(inout) :: A, B
  type(realization_type) :: realization
 ! integer, intent(inout) :: flag
  MatStructure flag

!   integer :: j, jn, jm1, jm2,jmu,mu, 
  integer :: ierr, i_upstream_revert
  integer :: nvar,neq,nr
  integer :: i1, i2, jng, i
  integer :: ithrm1, ithrm2
  integer :: ip1, ip2 
  integer :: p1,p2 
  real*8 :: dum1, dum2
  PetscViewer :: viewer

  PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  integer :: iicap,ii,jj,iiphas,iiphas1,iiphas2,iicap1,iicap2
  integer :: index_var_begin, index_var_end
! integer ibc_hencoeff
  real*8 :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho,dddt,dddp,fg,dfgdp,&
            dfgdt,eng,dhdt,dhdp,visc,dvdt,dvdp
! real*8 :: cond, gravity, acc, density_ave, 
  real*8 :: vv_darcy(realization%option%nphase), voldt, pvoldt
! real*8 :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
!           flux, fluxh, fluxv, difff, diffg, diffl,
  real*8 :: ff,dif(1:realization%option%nphase)
  real*8 :: tsrc1,qsrc1,csrc1
  real*8 :: dd1, dd2, dd, f1, f2
! real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
! real*8 :: por1, por2, den
  real*8 :: perm1, perm2
! real*8 :: qu_rate, p_vapor,sat_pressure_t
! real*8 :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
! real*8 :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! real*8 :: sat_pressure  ! Saturation pressure of water.
! real*8 :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu, sat_ave
  real*8 :: ra(1:realization%option%ndof,1:2*realization%option%ndof)  
! real*8 :: uhh, uconc, ukvr
! real*8 :: upweight,m1weight,m2weight,mbweight,mnweight
  real*8 :: delxbc(1:realization%option%ndof)
  real*8 :: blkmat11(1:realization%option%ndof,1:realization%option%ndof), &
            blkmat12(1:realization%option%ndof,1:realization%option%ndof),&
            blkmat21(1:realization%option%ndof,1:realization%option%ndof),&
            blkmat22(1:realization%option%ndof,1:realization%option%ndof)
  real*8 :: ResInc(1:realization%grid%nlmax, 1:realization%option%ndof, 1:realization%option%ndof),res(1:realization%option%ndof)  
  real*8 :: max_dev  

  integer :: local_id, ghosted_id
  integer :: local_id_up, local_id_dn
  integer :: ghosted_id_up, ghosted_id_dn
  integer ::  natural_id_up,natural_id_dn
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  logical :: enthalpy_flag
  integer :: iconn
  integer :: sum_connection  
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity, upweight
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option 
  type(field_type), pointer :: field  

!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  grid => realization%grid
  option => realization%option
  field => realization%field
  
  flag = SAME_NONZERO_PATTERN

! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
! call VecGetArrayF90(field%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p, ierr)

! print *,' In mph Jacobian ::  got pointers '
! ********************************************************************

! Accumulation terms

  ResInc=0.D0
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    
    voldt = volume_p(local_id) / option%dt
    pvoldt = porosity_loc_p(ghosted_id) * voldt

    iiphas=iphase_loc_p(ghosted_id)
 ! pressure equation    
    do nvar=1, option%ndof
   
      index_var_begin=(ghosted_id-1)*size_var_node+nvar*size_var_use+1
      index_var_end = index_var_begin -1 + size_var_use

      call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin : index_var_end),&
        porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(int(ithrm_loc_p(ghosted_id))),&
        option, Res,1,ierr)
      
      ResInc(local_id,:,nvar) = ResInc(local_id,:,nvar) + Res(:)
    enddo
  enddo

#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian1.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
! Source / Sink term
#if 0
  do nr = 1, option%nblksrc
      
    kk1 = option%k1src(nr) - option%nzs
    kk2 = option%k2src(nr) - option%nzs
    jj1 = option%j1src(nr) - option%nys
    jj2 = option%j2src(nr) - option%nys
    ii1 = option%i1src(nr) - option%nxs
    ii2 = option%i2src(nr) - option%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(option%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(option%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(option%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, option%ntimsrc
      if (option%timesrc(i,nr) == option%t) then
        tsrc1 = option%tempsrc(i,nr)
        qsrc1 = option%qsrc(i,nr)
        csrc1 = option%csrc(i,nr)
        goto 10
      else if (option%timesrc(i,nr) > option%t) then
        ff = option%timesrc(i,nr)-option%timesrc(i-1,nr)
        f1 = (option%t - option%timesrc(i-1,nr))/ff
        f2 = (option%timesrc(i,nr)-option%t)/ff
        tsrc1 = f1*option%tempsrc(i,nr) + f2*option%tempsrc(i-1,nr)
        qsrc1 = f1*option%qsrc(i,nr) + f2*option%qsrc(i-1,nr)
        csrc1 = f1*option%csrc(i,nr) + f2*option%csrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
   !print *,'pflow2ph : ', option%myrank,i,option%timesrc(i,nr), &
   !option%timesrc(i-1,nr),option%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
    
    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
            ng = grid%nL2G(n)
            
            do nvar=1,option%ndof      
              call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+2),&
              dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!             units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!             qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
              ResInc(n,option%jh2o,nvar)=  ResInc(n,option%jh2o,nvar) - qsrc1*option%dt
              ResInc(n,option%ndof,nvar)=  ResInc(n,option%ndof,nvar) - qsrc1*enth_src_h2o*option%dt

      
      
      !       print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*option%fmwh2o,qsrc1
            enddo
          enddo
        enddo
      enddo
    endif  
    
    if (csrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*option%nlx+(kk-1)*option%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*option%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            do nvar=1,option%ndof     
              rho = var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+4+option%nphase)*option%fmwco2 
              call co2_span_wagner(var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+2)*1.D-6,&
                tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)

         !    units: rho [kg/m^3]; csrc1 [kmol/s]

              enth_src_co2 = enth_src_co2 * option%fmwco2

              ResInc(n,option%jco2,nvar)=  ResInc(n,option%jco2,nvar) - csrc1*option%dt
              ResInc(n,option%ndof,nvar)=  ResInc(n,option%ndof,nvar) - csrc1*enth_src_co2*option%dt

          !   Res_AR(n,option%jco2)= Res_AR(n,option%jco2) - csrc1
    !         Res_AR(n,option%ndof)= Res_AR(n,option%ndof) - csrc1 * enth_src_co2
       !      r_p(s1) = r_p(s1) - csrc1

!             print *,'pflow2ph_co2: ',nr,n,ng,tsrc1,rho,option%fmwco2,csrc1
            enddo
          enddo
        enddo
      enddo
    endif
  enddo  
#endif
  
  ! print *,' Mph Jaco Finished source terms'
  
! Contribution from BC
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection

    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if(ghosted_id<=0)then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      end if
  
      p1 = 1 + (local_id-1) * option%ndof
      
      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_object%dist(3,iconn))
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = cur_connection_object%dist(3,iconn) * &
                         cur_connection_object%dist(0,iconn)
       
      delxbc=0.D0
      select case(boundary_condition%condition%itype(1))
        case(1)
          delxbc=0.D0
        case(2)
          ! solve for pb from Darcy's law given qb /= 0
         field%xxbc(:,sum_connection)=xx_loc_p((ghosted_id-1)*option%ndof+1: ghosted_id*option%ndof)
          field%iphasebc(sum_connection) = int(iphase_loc_p(ghosted_id))
             delxbc=option%delx(1:option%ndof,ghosted_id)
        case(3) 
          ! field%xxbc(1,iconn)=option%pressurebc(2,ibc)
          field%xxbc(2:option%ndof,sum_connection)= xx_loc_p((ghosted_id-1)*option%ndof+2: ghosted_id*option%ndof)
          field%iphasebc(sum_connection)=int(iphase_loc_p(ghosted_id))
          delxbc(1)=0.D0
          delxbc(2:option%ndof)=option%delx(2:option%ndof,ghosted_id)
        case(4)
          field%xxbc(1,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+1)
          field%xxbc(3:option%ndof,sum_connection) = xx_loc_p((ghosted_id-1)*option%ndof+3: ghosted_id*option%ndof)    
          delxbc(1)=option%delx(1,ghosted_id)
          delxbc(3:option%ndof) = option%delx(3:option%ndof,ghosted_id) 
          field%iphasebc(sum_connection)=int(iphase_loc_p(ghosted_id))

      end select

! print *,'2ph bc',option%myrank,nc,m,ng,ibc,option%ibndtyp(ibc),option%pressurebc(:,ibc), &
! option%tempbc(ibc),option%sgbc(ibc),option%concbc(ibc),field%velocitybc(:,ibc)

!   if(option%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   option%ibndphaseRate(ibc) 
   !   option%ibndconc(ibc)    ! 
   !   option%tempbc(ibc)      !1 elements 
   !   option%pressurebc(ibc)  !nphase elements
!      endif
   
      iicap = int(icap_loc_p(ghosted_id))     
       
!      print *,'pflow_2pha_bc: ',option%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      option%nconnbc,option%ibndtyp(ibc),option%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,option%ideriv,option%ibndtyp(ibc),option%density_bc,&
!   option%pressurebc(2,ibc),option%tempbc(ibc),option%concbc(ibc),option%sgbc(ibc)
        !*****************
      dif(1)= option%difaq
      dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    !*******************************************

  ! here should pay attention to BC type !!!
      call pri_var_trans_mph_ninc(field%xxbc(:,sum_connection),field%iphasebc(sum_connection),&
                              option%scale,option%nphase,option%nspec,iicap,dif,&
                              field%varbc(1:size_var_use),option%itable, &
                              option%m_nacl,ierr, dum1, dum2)
  
      call pri_var_trans_mph_winc(field%xxbc(:,sum_connection),delxbc,&
                              field%iphasebc(sum_connection),option%scale,option%nphase, &
                              option%nspec,iicap,dif(1:option%nphase), &
                              field%varbc(size_var_use+1: &
                                         (option%ndof+1)*size_var_use), &
                              option%itable,option%m_nacl,ierr)
            
!    print *,' Mph Jaco BC terms: finish increment'
      do nvar=1,option%ndof
   
        call MPHASERes_FLBCCont(sum_connection,boundary_condition%condition%itype(1), &
                                cur_connection_object%area(iconn), &
                                field%varbc(nvar*size_var_use+1: &
                                           (nvar+1)*size_var_use), &
                                var_loc_p((ghosted_id-1)*size_var_node+ &
                                             nvar*size_var_use+1: &
                                          (ghosted_id-1)*size_var_node+ &
                                             nvar*size_var_use+size_var_use), &
                                porosity_loc_p(ghosted_id),tor_loc_p(ghosted_id), &
                                option%sir(1:option%nphase,iicap), &
                                cur_connection_object%dist(0,iconn),&
                                perm1,D2,distance_gravity,option,field,vv_darcy,Res)

    
        ResInc(local_id,1:option%ndof,nvar) = ResInc(local_id,1:option%ndof,nvar) - Res(1:option%ndof)
      enddo

    enddo
    boundary_condition => boundary_condition%next
  enddo

  do local_id= 1, grid%nlmax
    ra=0.D0
    
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    
    natural_id_up= grid%nG2N(ghosted_id)
   ! Remember, the matrix index starts from (0,0)
    p1 = (ghosted_id-1)*option%ndof ! = 1 + (ng-1)*option%ndof-1
   
    max_dev=0.D0
    do neq=1, option%ndof
      do nvar=1, option%ndof
        ra(neq,nvar)=ResInc(local_id,neq,nvar)/option%delx(nvar,ghosted_id) - &
                     ResOld_AR(local_id,neq)/option%delx(nvar,ghosted_id)
          if(max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
     
      enddo      
    enddo
    if(option%use_isoth==PETSC_TRUE)then
      ra(:,2)=0.D0
      ra(3,1:option%ndof)=0.D0
      ra(3,2)=1.D0
    endif     

   
  ! if(max_dev<1D-5)then
  !  print *,'Mph Jaco max dev = ', max_dev
 !  endif
  
    select case(option%idt_switch)
      case(1) 
        ra(1:option%ndof,1:option%ndof) =ra(1:option%ndof,1:option%ndof) /option%dt
      case(-1)
        if(option%dt>1) ra(1:option%ndof,1:option%ndof) =ra(1:option%ndof,1:option%ndof) /option%dt
    end select
      

    if (option%iblkfmt == 0) then
      p1=(natural_id_up)*option%ndof
     
      do ii=0,option%ndof-1
        do jj=0,option%ndof-1
          call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1)/ volume_p(local_id),ADD_VALUES,ierr)
        enddo
       enddo
    else
      !ra(1:option%ndof,1:option%ndof) =ra(1:option%ndof,1:option%ndof) /option%dt
      blkmat11=ra(1:option%ndof,1:option%ndof)
    
      if(volume_p(local_id)>1.D0 ) blkmat11=blkmat11 / volume_p(local_id)
   
     ! if(n==1) print *,  blkmat11, volume_p(n), ra
      call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_up,blkmat11,ADD_VALUES,ierr)
    endif
         
  enddo
!   print *,' Mph Jaco Finished one node terms'
! -----------------------------contribution from transport----------------------
#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian2.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif


 !print *,'phase cond: ',iphase_loc_p
  ResInc=0.D0
 
  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_object%id_up(iconn)
      ghosted_id_dn = cur_connection_object%id_dn(iconn)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or. field%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      natural_id_up = grid%nG2N(ghosted_id_up)
      natural_id_dn = grid%nG2N(ghosted_id_dn)
      p2 =  (ghosted_id_dn-1)*option%ndof
   
      fraction_upwind = cur_connection_object%dist(-1,iconn)
      distance = cur_connection_object%dist(0,iconn)
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = cur_connection_object%dist(3,iconn)*distance
      dd1 = distance*fraction_upwind
      dd2 = distance-dd1 ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd2/(dd1+dd2)
    
      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_object%dist(3,iconn))

      perm2 = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_object%dist(3,iconn))
    
      iiphas1 = iphase_loc_p(ghosted_id_up)
      iiphas2 = iphase_loc_p(ghosted_id_dn)

      ithrm1 = ithrm_loc_p(ghosted_id_up)
      ithrm2 = ithrm_loc_p(ghosted_id_dn)
      D1 = option%ckwet(ithrm1)
      D2 = option%ckwet(ithrm2)
    
      iicap1 = int(icap_loc_p(ghosted_id_up))
      iicap2 = int(icap_loc_p(ghosted_id_dn))
  
  ! do neq = 1, option%ndof
      do nvar = 1, option%ndof
    
        call MPHASERes_FLCont(sum_connection,cur_connection_object%area(iconn), &
                              var_loc_p((ghosted_id_up-1)*size_var_node+ &
                                          nvar*size_var_use+1: &
                                        (ghosted_id_up-1)*size_var_node+ &
                                          nvar*size_var_use+size_var_use), &
                              porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                              option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                              var_loc_p((ghosted_id_dn-1)*size_var_node+1: &
                                        (ghosted_id_dn-1)*size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                              option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                              distance_gravity,upweight, &
                              option, vv_darcy,Res)

        ra(:,nvar)= Res(:)/option%delx(nvar,ghosted_id_up)-ResOld_FL(sum_connection,:)/option%delx(nvar,ghosted_id_up)

!     if(vv_darcy(1)>0.D0 .and. option%iupstream(sum_connection,1) == -1) i_upstream_revert =1 
!     if(vv_darcy(1)<0.D0 .and. option%iupstream(sum_connection,1) == 1)  i_upstream_revert =1
!     if(vv_darcy(2)>0.D0 .and. option%iupstream(sum_connection,2) == -1) i_upstream_revert =2
!     if(vv_darcy(2)<0.D0 .and. option%iupstream(sum_connection,2) == 1 ) i_upstream_revert =2
  
       
        call MPHASERes_FLCont(sum_connection,cur_connection_object%area(iconn), &
                              var_loc_p((ghosted_id_up-1)*size_var_node+1: &
                                        (ghosted_id_up-1)*size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                              option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                              var_loc_p((ghosted_id_dn-1)*size_var_node+ &
                                          nvar*size_var_use+1: &
                                        (ghosted_id_dn-1)*size_var_node+ &
                                          nvar*size_var_use+size_var_use), &
                              porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                              option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                              distance_gravity,upweight, &
                              option, vv_darcy,Res)
 
        ra(:,nvar+option%ndof)= Res(:)/option%delx(nvar,ghosted_id_dn)-ResOld_FL(sum_connection,:)/option%delx(nvar,ghosted_id_dn)

     
!     if(vv_darcy(1)>0.D0 .and. option%iupstream(sum_connection,1) == -1) i_upstream_revert =3 
!     if(vv_darcy(1)<0.D0 .and. option%iupstream(sum_connection,1) == 1)  i_upstream_revert =3
!     if(vv_darcy(2)>0.D0 .and. option%iupstream(sum_connection,2) == -1) i_upstream_revert =4
!     if(vv_darcy(2)<0.D0 .and. option%iupstream(sum_connection,2) == 1 ) i_upstream_revert =4
   
      enddo
  
   !   print *,' Mph Jaco Finished NC terms'
  
  ! enddo
   
      if (option%use_isoth==PETSC_TRUE) then
        ra(3,1:2*option%ndof)=0.D0
        ra(:,2)=0.D0
        ra(:,2+option%ndof)=0.D0
      endif   
 

      if (option%iblkfmt == 1) then
        blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
      endif
    
      p1=(natural_id_up)*option%ndof
      p2=(natural_id_dn)*option%ndof
 
      select case(option%idt_switch)
        case(1)
          ra =ra / option%dt
        case(-1)  
          if (option%dt>1) ra =ra / option%dt
      end select  
       
      do ii=0,option%ndof-1
        do jj=0,option%ndof-1
          if (local_id_up>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1)/volume_p(local_id_up),ADD_VALUES,ierr)
            else
              blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
            endif
          endif
          if (local_id_dn>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p2+ii,p1+jj,-ra(ii+1,jj+1)/volume_p(local_id_dn),ADD_VALUES,ierr)
            else
              blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
            endif
          endif
        enddo
   
        do jj=option%ndof,2*option%ndof-1
          if (local_id_up>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p1+ii,p2+jj-option%ndof,ra(ii+1,jj+1)/volume_p(local_id_up),ADD_VALUES,ierr)
            else
              blkmat12(ii+1,jj-option%ndof+1) = blkmat12(ii+1,jj-option%ndof+1) + ra(ii+1,jj+1)
            endif
          endif
          if (local_id_dn>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p2+ii,p2+jj-option%ndof,-ra(ii+1,jj+1)/volume_p(local_id_dn),ADD_VALUES,ierr)
            else
              blkmat22(ii+1,jj-option%ndof+1) =  blkmat22(ii+1,jj-option%ndof+1) - ra(ii+1,jj+1)
            endif
          endif
        enddo
      enddo
  
      if (option%iblkfmt /= 0) then
        if (volume_p(local_id_up)>1.D0) then
          blkmat11=blkmat11/volume_p(local_id_up); blkmat12=blkmat12/volume_p(local_id_up)
        endif
        if (volume_p(local_id_dn)>1.D0) then
          blkmat21=blkmat21/volume_p(local_id_dn); blkmat22=blkmat22/volume_p(local_id_dn)
        endif 
   
       !  if(dabs(volume_p(n1)-3.D0)>1D-5 .and. n1>0) print *, n1,  volume_p(n1)
   !   if(dabs(volume_p(n2)-3.D0)>1D-5 .and. n2>0) print *, n2,  volume_p(n2)
        if(local_id_up>0)call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_up,blkmat11,ADD_VALUES,ierr)
        if(local_id_dn>0)call MatSetValuesBlocked(A,1,natural_id_dn,1,natural_id_dn,blkmat22,ADD_VALUES,ierr)
        if(local_id_up>0)call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_dn,blkmat12,ADD_VALUES,ierr)
        if(local_id_dn>0)call MatSetValuesBlocked(A,1,natural_id_dn,1,natural_id_up,blkmat21,ADD_VALUES,ierr)
      endif
!print *,'accum r',ra(1:5,1:8)   
 !print *,'devq:',nc,q,dphi,devq(3,:)
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
 

  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%var_loc, var_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  f1 = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,f1,ierr) 
  do i=1, n_zero_rows
    n = mod(zero_rows_local(i),option%ndof)
    p1 = zero_rows_local_ghosted(i)
    if (n == 0) then
      p2 = p1-1
    if (n == option%ndof-1) then
      p2 = p1+1
    else
      p2 = p1
    endif
    call MatSetValuesLocal(A,1,p1,1,p2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  f1 = 1.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,f1,ierr) 
#endif

#ifdef DEBUG_GEH    
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine MPHASEJacobian



subroutine pflow_mphase_initaccum(realization)
 
  use translator_mph_module, only : pri_var_trans_mph_ninc
  
  use Realization_module
  use Grid_module
  use Field_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization 
  
  integer :: ierr
  integer :: i, index_var_begin,index_var_end
  integer :: p1
  integer :: iicap, iiphase
  integer :: local_id, ghosted_id  

  PetscScalar, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_loc_p(:),&
                          var_loc_p(:), icap_loc_p(:),iphase_loc_p(:),ithrm_loc_p(:)
  
 !integer, pointer ::iphase_loc_p(:)
  
! real*8 :: sat_pressure
  real*8 :: pvol, satw  ! Saturation pressure of water.
  real*8 :: dif(1:realization%option%nphase),res(1:realization%option%ndof)
 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%accum, accum_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  !print *,'mphase initaccum  Gotten pointers'
 
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    iicap=int(icap_loc_p(ghosted_id))
    iiphase = int(iphase_loc_p(ghosted_id))
    dif(1)= option%difaq
    dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))

    call pri_var_trans_mph_ninc(yy_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase,&
                                option%scale,option%nphase,option%nspec, iicap, dif,&
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),option%itable,&
                                option%m_nacl,ierr, satw, pvol)

  enddo

  call VecRestoreArrayF90(field%var_loc, var_loc_p,ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p,ierr)

!---------------------------------------------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    p1 = 1 + (local_id-1)*option%ndof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    i = ithrm_loc_p(ghosted_id)
    
    call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin: index_var_end),&
                          porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(i), option, &
                          Res, 0,ierr)
 

    accum_p(p1:p1+option%ndof-1)=Res(:) 

   !print *, 'init m accum ', n,  Res 

! print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*option%nspec+(n-1)*option%nphase*option%nspec),kvr_p(jn),kvr_p(jn+1)

  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%accum, accum_p, ierr)
  call VecRestoreArrayF90(field%var_loc, var_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine pflow_mphase_initaccum


subroutine pflow_update_mphase(realization)
  
  use translator_mph_module, only : pri_var_trans_mph_ninc, translator_mph_get_output, translator_mphase_massbal
  
  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module  

  implicit none

  type(realization_type) :: realization 
    
  integer :: dof_offset
    integer :: iithrm
    integer :: ierr,iicap,iiphase
    PetscScalar, pointer :: xx_p(:),icap_loc_p(:),ithrm_loc_p(:),iphase_loc_p(:), var_loc_p(:)
    real*8 :: dif(1:realization%option%nphase), dum1, dum2           
  integer :: local_id, ghosted_id     

  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  integer :: sum_connection  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => realization%grid
  option => realization%option
  field => realization%field
        
  ! if (option%rk > 0.d0) call Rock_Change(grid)
  ! call Translator_MPhase_Switching(field%xx,grid,1,ierr)
  ! print *,'MPhase_Update done'
 
   !if(ichange ==1)then
    call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
    call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
    call VecGetArrayF90(field%var_loc,var_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)    
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
    
    iicap = icap_loc_p(ghosted_id)
    iiphase = iphase_loc_p(ghosted_id)
    dof_offset=(local_id-1)*option%ndof
    if(xx_p(dof_offset+3)<0.D0) xx_p(dof_offset+3)=zerocut

    ! *****************
      dif(1)= option%difaq
      dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    ! *******************************************
      call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase,&
      option%scale,option%nphase,option%nspec, iicap, dif,&
      var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),&
      option%itable,option%m_nacl,ierr, dum1, dum2)

    enddo
#if 0
  !geh added for transient boundary conditions  
  if (associated(field%imat) .and. option%iread_geom < 0) then

    boundary_condition => realization%boundary_conditions%first
    sum_connection = 0
    do 
      if (.not.associated(boundary_condition)) exit
    
      cur_connection_object => boundary_condition%connection

      do iconn = 1, cur_connection_object%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_object%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        if (associated(field%imat)) then
          if (field%imat(ghosted_id) <= 0) cycle
        endif
       
        if (local_id<0) then
          print *, "Wrong boundary node index... STOP!!!"
          stop
        endif

        if (boundary_condition%condition%itype(1)==1 .or. &
            boundary_condition%condition%itype(1)==3) then
          iicap=int(icap_loc_p(ghosted_id))
          iithrm=int(ithrm_loc_p(ghosted_id)) 
          dif(1)= option%difaq
      
          if(field%iphasebc(nc) ==3)then
            sw= field%xxbc(1,sum_connection)
            call pflow_pckr_richards_fw(iicap ,sw,pc,kr)    
            field%xxbc(1,sum_connection) =  option%pref - pc(1)
          endif
      
          call pri_var_trans_Richards_ninc(field%xxbc(:,sum_connection),field%iphasebc(sum_connection), &
                                           option%scale,option%nphase,option%nspec, &
                                           iicap,dif, &
                                           field%varbc(1:size_var_use),option%itable,ierr, &
                                           option%pref)
        
          if (translator_check_cond_Richards(field%iphasebc(sum_connection), &
                                             field%varbc(1:size_var_use), &
                                             option%nphase,option%nspec) /=1) then
            print *," Wrong bounday node init...  STOP!!!", field%xxbc(:,sum_connection)
      
            print *,field%varbc
            stop    
          endif 
        endif

        if (boundary_condition%condition%itype(1)==2) then
          yybc(2:option%ndof,nc)= field%xxbc(2:option%ndof,nc)
          vel_bc(1,nc) = field%velocitybc(1,nc)
        endif 
      
      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif
#endif
   
    call VecRestoreArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
    call VecRestoreArrayF90(field%var_loc,var_loc_p,ierr)
     
    if(option%nphase>1) call translator_mphase_massbal(realization)
   ! endif 

    call VecCopy(field%xx, field%yy, ierr)   
    call VecCopy(field%iphas_loc, field%iphas_old_loc, ierr)   
     
    call  pflow_mphase_initaccum(realization)
      !print *,'pflow_mphase_initaccum done'
    call translator_mph_get_output(realization)
 !  print *,'translator_get_output done'
  ! the output variables should be put into option%pressure, temp,xmol,sat...
  ! otherwise need to rewrite the pflow_output

end subroutine pflow_update_mphase


subroutine pflow_mphase_initadj(realization)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  ! geh - will not compile without the 'only:' statement
  use translator_mph_module, only : pri_var_trans_mph_ninc, translator_check_phase_cond

  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module
  implicit none
 
  type(realization_type) :: realization 

  integer :: ierr
  integer :: num_connection  
  integer :: jn
  integer :: iicap
  integer :: iiphase,iithrm
  integer :: local_id, ghosted_id 

  PetscScalar, pointer :: xx_p(:),var_loc_p(:)
  PetscScalar, pointer ::iphase_loc_p(:), ithrm_loc_p(:),icap_loc_p(:)
  
  real*8  dif(realization%option%nphase), dum1, dum2
  
! real*8 :: temp1
! real*8, parameter :: Rg=8.31415D0


  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  integer :: sum_connection
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => realization%grid
  option => realization%option 
  field => realization%field 


  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%xx, xx_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p, ierr)
! print *,'initadj gotten pointers' 


 do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
     
    jn = 1 + (local_id-1)*option%nphase
    iicap=int(icap_loc_p(ghosted_id))
        
    iiphase = iphase_loc_p(ghosted_id)
    !*****************
    dif(1)= option%difaq
    dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    !******************************************* 
    call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase,&
      option%scale,option%nphase,option%nspec, iicap,  dif,&
      var_loc_p((ghosted_id-1)*size_var_node+1: (ghosted_id-1)*size_var_node+size_var_use), &
      option%itable,option%m_nacl,ierr, dum1, dum2)
   
   !print *, xx_p((n-1)*option%ndof+1:n*option%ndof)
    if(translator_check_phase_cond(iiphase, &
      var_loc_p((ghosted_id-1)*size_var_node+1: (ghosted_id-1)*size_var_node+size_var_use),&
      option%nphase,option%nspec) /= 1 ) then
      print *," Wrong internal node init...  STOP!!!"
      stop    
    endif 
  enddo


  boundary_condition => realization%boundary_conditions%first
  num_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit    
    num_connection = num_connection + boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo

! implemented in richards, but not here????
!  allocate(yybc(option%ndof,num_connection))
!  allocate(vel_bc(option%nphase,num_connection))
!  yybc =field%xxbc
!  vel_bc = field%velocitybc

  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0  
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection
    
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
  
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif
       

      if(local_id<0)then
         print *, "Wrong boundary node index... STOP!!!"
         stop
      end if

      if (boundary_condition%condition%itype(1)==1) then
!          boundary_condition%condition%itype(1)==3) then      
        iicap=int(icap_loc_p(ghosted_id))
        iithrm=int(ithrm_loc_p(ghosted_id)) 
        dif(1)= option%difaq
        dif(2)= option%cdiff(iithrm)

        call pri_var_trans_mph_ninc(field%xxbc(:,sum_connection),field%iphasebc(sum_connection),&
                                    option%scale,option%nphase,option%nspec,iicap, &
                                    dif,field%varbc(1:size_var_use), &
                                    option%itable,option%m_nacl,ierr, dum1, dum2)
      
        if (translator_check_phase_cond(field%iphasebc(sum_connection), &
                                        field%varbc(1:size_var_use), &
                                        option%nphase,option%nspec) &
            /=1) then
          print *," Wrong bounday node init...  STOP!!!", field%xxbc(:,iconn)
      
          print *,field%varbc
          stop    
        endif 
      endif

#if 0  
! not implemented in co2???
      if (boundary_condition%condition%itype(1)==2) then
    
        yybc(2:option%ndof,nc)= field%xxbc(2:option%ndof,nc)
        vel_bc(1,nc) = field%velocitybc(1,nc)
      endif 
#endif      
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%xx, xx_p, ierr)
  call VecRestoreArrayF90(field%var_loc, var_loc_p, ierr)
  !print *,kgjkdf
  
  !call VecCopy(field%iphas,field%iphas_old,ierr)

end subroutine pflow_mphase_initadj



subroutine createmphaseZeroArray(realization)

  use Realization_module
  use Grid_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization
  integer :: ncount, idof
  integer :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
    
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  n_zero_rows = 0

  if (associated(field%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (field%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%ndof
      else
#ifdef ISOTHERMAL
        n_zero_rows = n_zero_rows + 1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    n_zero_rows = n_zero_rows + grid%nlmax
#endif
  endif

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))
  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  if (associated(field%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (field%imat(ghosted_id) <= 0) then
        do idof = 1, option%ndof
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%ndof+idof
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%ndof+idof-1
        enddo
      else
#ifdef ISOTHERMAL
        ncount = ncount + 1
        zero_rows_local(ncount) = local_id*option%ndof
        zero_rows_local_ghosted(ncount) = ghosted_id*option%ndof-1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ncount = ncount + 1
      zero_rows_local(ncount) = local_id*option%ndof
      zero_rows_local_ghosted(ncount) = ghosted_id*option%ndof-1
    enddo
#endif
  endif

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine createmphaseZeroArray

end module MPHASE_module
