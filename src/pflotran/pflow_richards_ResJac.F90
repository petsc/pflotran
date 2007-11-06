! introduced grid variables: e_total :: 1 dof
!Richards translator_module                           c_total :: option%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (option%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

               
module Richards_module

  private 
#include "include/finclude/petsc.h"
!#include "include/petscf90.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
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
  real*8, parameter :: formeps   = 100.D0
  real*8, parameter :: eps       = 1.D-8
  real*8, parameter :: floweps   = 1.D-24
  real*8, parameter :: satcuteps = 1.D-8
  real*8, parameter :: dfac = 1.D-8

  integer,save :: size_var_use 
  integer,save :: size_var_node
  real*8, allocatable,save :: Resold_AR(:,:), Resold_FL(:,:)
  real*8, pointer, save :: yybc(:,:), vel_bc(:,:)
! Contributions to residual from accumlation/source/Reaction, flux(include diffusion)
  
  
  
   

  public RichardsResidual, RichardsJacobian, pflow_Richards_initaccum, &
         pflow_update_Richards,pflow_Richards_initadj, pflow_Richards_timecut,&
         pflow_Richards_setupini, Richards_Update, Richards_Update_Reason

  public :: createRichardsZeroArray
  integer, save :: n_zero_rows = 0
  integer, pointer, save :: zero_rows_local(:)  ! 1-based indexing
  integer, pointer, save :: zero_rows_local_ghosted(:) ! 0-based indexing

contains


subroutine pflow_Richards_timecut(solution)
 
  use Solution_module
  use Option_module
  use Grid_module
 
  implicit none
  
  type(solution_type) :: solution
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  
  PetscScalar, pointer :: xx_p(:),yy_p(:)!,var_p(:),iphase_p(:)
  integer :: n,n0,re,ierr
  !integer re0, ierr, index, iiphaRichards
  !real*8, pointer :: sat(:),xmol(:)

  grid => solution%grid
  option => solution%option
 
  call VecGetArrayF90(option%xx, xx_p, ierr)
  call VecGetArrayF90(option%yy, yy_p, ierr)
 ! call VecGetArrayF90(option%var, var_p, ierr); 
 ! call VecGetArrayF90(option%iphas, iphase_p, ierr); 

  do n=1, grid%nlmax
    n0=(n-1)*option%ndof
    do re = 1, option%ndof
   !xx_p(n0+re)= 0.5D0 * xx_p(n0+re) +.5D0 *yy_p(n0+re)
      xx_p(n0+re)= yy_p(n0+re)
    enddo
  enddo 
  call VecRestoreArrayF90(option%xx, xx_p, ierr) 
  call VecRestoreArrayF90(option%yy, yy_p, ierr)
  
  !call VecCopy(option%xx,option%yy,ierr)
  !call pflow_Richards_initaccum(grid)
 
end subroutine pflow_Richards_timecut
  

subroutine pflow_Richards_setupini(solution)

  use Solution_module
  use Option_module
  use Grid_module
  use Region_module
  use Structured_Grid_module
  use Coupler_module
  use Condition_module
 
  implicit none
  
  type(solution_type) :: solution
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: initial_condition

  PetscScalar, pointer :: xx_p(:), iphase_p(:)
  integer local_id, ibegin, iend, icell, ierr
  
  grid => solution%grid
  option => solution%option
  
  size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
  size_var_node = (option%ndof + 1) * size_var_use
  
  allocate(Resold_AR(grid%nlmax,option%ndof))
  allocate(Resold_FL(grid%internal_connection_list%first%num_connections, &
                     option%ndof))
  allocate(option%delx(option%ndof,grid%ngmax))
  option%delx=0.D0
   
  call VecGetArrayF90(option%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%iphas, iphase_p,ierr)
  
  initial_condition => solution%initial_conditions%first
  
  do
  
    if (.not.associated(initial_condition)) exit
    
    do icell=1,initial_condition%region%num_cells
      local_id = initial_condition%region%cell_ids(icell)
      iend = local_id*option%ndof
      ibegin = iend-option%ndof+1
      xx_p(ibegin:iend) = &
        initial_condition%flow_condition%cur_value(1:option%ndof)
      iphase_p(local_id)=initial_condition%flow_condition%iphase
    enddo
  
    initial_condition => initial_condition%next
  
  enddo
  
  call VecRestoreArrayF90(option%xx, xx_p, ierr)
  call VecRestoreArrayF90(option%iphas, iphase_p,ierr)

end  subroutine pflow_Richards_setupini
  

subroutine Richards_Update_Reason(reason,solution)

  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none
 
  type(solution_type) :: solution
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  integer, intent(out):: reason
  PetscScalar, pointer :: xx_p(:),var_p(:),iphase_p(:), yy_p(:) !,r_p(:)
  integer :: n,n0,re
  integer re0, ierr, iipha
! integer :: index
  
  grid => solution%grid
  option => solution%option
  
! real*8, pointer :: sat(:),xmol(:)
! real*8 rmax(option%ndof)

! why this barrier
  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  re = 1
 ! call SNESComputeFunction(option%snes,option%xx,grid%r,ierr)
 ! do n=1,option%ndof
 !  call VecStrideNorm(grid%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if (re>0) then
    call VecGetArrayF90(option%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(option%yy, yy_p, ierr)
    call VecGetArrayF90(option%var, var_p, ierr); 
    call VecGetArrayF90(option%iphas, iphase_p, ierr); 
  
    do n = 1,grid%nlmax

      !geh - Ignore inactive cells with inactive materials
      if (associated(option%imat)) then
        if (option%imat(grid%nL2G(n)) <= 0) cycle
      endif

      n0=(n-1)* option%ndof
      !index=(n-1)*size_var_node
      !sat=>var_p(index+2+1:index+2+option%nphase)
      !den=>var_p(index+2+option%nphase+1:index+2+2*option%nphase)
    !xmol=>var_p(index+2+7*option%nphase+1:index+2+7*option%nphase + option%nphase*option%nspec)    
      iipha=int(iphase_p(n))
     !if(n==3583 .or. n==3587)
   !print *, 'update reson', grid%nlmax, n, iipha, xx_p(n0+1:n0+3)
   !if(xmol(4)>1.0) re=0; goto 1
   !if(xmol(4)<.0) re=0; goto 1
   !if(sat(2) < .0) re=0;goto 1
   ! if(sat(2) > 1.) re=0;goto 1
      select case(iipha)
        case (1)
          if (xx_p(n0 + 1) < option%pref) then
            re = 0
            exit
          endif
        case (3)
          if (xx_p(n0 + 1) > option%pref) then
            re = 0
            exit
          endif
      end select  
    enddo
  
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
   ! print *, 'update reason: ',grid%myrank,grid%nlmax,n,re

  !   call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
   !print *,' update reason ba MPI', ierr
    if (re<=0) print *,'Sat or Con out of Region at: ',n,iipha,xx_p(n0+1:n0+2)
    call VecRestoreArrayF90(option%xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(option%yy, yy_p, ierr)
    call VecRestoreArrayF90(option%var, var_p, ierr) 
    call VecRestoreArrayF90(option%iphas, iphase_p, ierr) 
  endif
 ! print *,' update reason', grid%myrank, re,n,grid%nlmax

  
  if (option%commsize >1) then
    call MPI_ALLREDUCE(re, re0,1, MPI_INTEGER,MPI_SUM, &
                       PETSC_COMM_WORLD,ierr)
  !print *,' update reason re'
    !call MPI_BCAST(re0,1, MPI_INTEGER, 0,PETSC_COMM_WORLD,ierr)
  !print *,' update reason ca'
    if (re0<option%commsize) re = 0
  endif
  reason = re
  
  if (reason<=0) print *,'Sat or Con out of Region'
  
end subroutine Richards_Update_Reason

!=======================================================================================
 


subroutine RichardsRes_ARCont(node_no, var_node,por,vol,rock_dencpr, option, Res_AR,ireac,ierr)

  use Option_module
  
  implicit none

  type(option_type) :: option
  integer node_no
  integer, optional:: ireac,ierr
  real*8, target:: var_node(1:size_var_use)
  real*8 Res_AR(1:option%ndof) 
  real*8 vol,por,rock_dencpr
     
  real*8, pointer :: temp, pre_ref   ! 1 dof
  real*8, pointer :: sat(:), density(:), amw(:), h(:), u(:), pc(:), kvr(:)         ! nphase dof
  real*8, pointer :: xmol(:), diff(:)            ! nphase*nspec
  
  integer :: ibase, m,np, iireac=1
  real*8 pvol,mol(option%nspec),eng
  
  if (present(ireac)) iireac=ireac
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
  ibase=ibase+option%nphase*option%nspec; diff=>var_node(ibase:ibase+option%nphase*option%nspec-1)

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
  if (iireac>0) then
!H2O
    mol(1)= mol(1) - option%dt * option%rtot(node_no,1)
!CO2
!    mol(2)= mol(2) - option%dt * option%rtot(node_no,2)
!should include related energy change here
  endif
  Res_AR(1:option%ndof-1)=mol(:)
  Res_AR(option%ndof)=eng
  nullify(temp, pre_ref, sat, density, amw, h,u, pc,kvr,xmol,diff)       

end subroutine  RichardsRes_ARCont


#ifdef OVERHAUL
#ifdef OVERHAUL_V2
subroutine RichardsRes_FLCont(nconn_no,area,var_node1,por1,tor1,sir1,dd1,perm1, &
                            Dk1,var_node2,por2,tor2,sir2,dd2,perm2,Dk2, &
                            dist_gravity,upweight, &
                            option,vv_darcy,Res_FL)
  use Option_module                              
#else
subroutine RichardsRes_FLCont(nconn_no,area,var_node1,por1,tor1,sir1,dd1,perm1, &
                            Dk1,var_node2,por2,tor2,sir2,dd2,perm2,Dk2,grid, &
                            vv_darcy,connection_object,Res_FL)
#endif
  use Connection_module
#else
subroutine RichardsRes_FLCont(nconn_no,area,var_node1,por1,tor1,sir1,dd1,perm1, &
                            Dk1,var_node2,por2,tor2,sir2,dd2,perm2,Dk2,option, &
                            vv_darcy,Res_FL)
#endif

  implicit none
  
  integer nconn_no
  type(option_type) :: option
  real*8 sir1(1:option%nphase),sir2(1:option%nphase)
  real*8, target:: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8, target:: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8 por1,por2,tor1,tor2,perm1,perm2,Dk1,Dk2,dd1,dd2
  real*8 vv_darcy(option%nphase),area
  real*8 Res_FL(1:option%ndof) 
  
  real*8 :: dist_gravity  ! distance along gravity vector
     
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  real*8, pointer :: xmol1(:), diff1(:)            ! 
  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  real*8, pointer :: xmol2(:), diff2(:)    
  
  integer ibase, m,np, ind
  real*8  fluxm(option%nspec),fluxe, v_darcy,q
  real*8 uh,uxmol(1:option%nspec), ukvr,difff,diffdp, DK,Dq
  real*8 upweight,density_ave,cond, gravity, dphi

#ifndef OVERHAUL_V2
#ifdef OVERHAUL  
  type(connection_type) :: connection_object
#endif
#endif
  
!  m1=grid%nd1(nc); n1 = grid%nG2L(m1) ! = zero for ghost nodes 
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
  ibase=ibase+option%nphase; 
                         xmol1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)
                         xmol2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)
  ibase=ibase+option%nphase*option%nspec;
                         diff1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)
                         diff2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)

  !print *,' FLcont got pointers' ,var_node1,var_node2,sir1,sir2
  !print *,' tmp=',temp1,temp2
  !print *,'diff=',diff1,diff2
   
  Dq = (perm1 * perm2)/(dd1*perm2 + dd2*perm1)
  diffdp = (por1 *tor1 * por2*tor2) / (dd2*por1*tor1 + dd1*por2*tor2)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy = 0.D0  
  
  do np=1, option%nphase

! Flow term
    if ((sat1(np) > sir1(np)) .or. (sat2(np) > sir2(np))) then
#ifndef OVERHAUL_V2    
      upweight=dd1/(dd1+dd2)
#endif      
      if (sat1(np) <eps) then 
        upweight=0.d0
      else if (sat2(np) <eps) then 
        upweight=1.d0
      endif
      density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  

#ifndef OVERHAUL    
      gravity = (upweight*density1(np)*amw1(np) + &
                (1.D0-upweight)*density2(np)*amw2(np)) &
                * option%gravity * grid%delz(nconn_no) * grid%grav_ang(nconn_no)
#else
#ifdef OVERHAUL_V2
      gravity = (upweight*density1(np)*amw1(np) + &
                (1.D0-upweight)*density2(np)*amw2(np)) &
                * option%gravity * dist_gravity
#else
      gravity = (upweight*density1(np)*amw1(np) + &
                (1.D0-upweight)*density2(np)*amw2(np)) &
                * option%gravity * connection_object%dist(3,nconn_no) * &
                connection_object%dist(0,nconn_no) 
#endif                
#endif

      dphi = pre_ref1 - pre_ref2  + gravity
!    print *,'FLcont  dp',dphi
  ! note uxmol only contains one phase xmol
      if (dphi>=0.D0) then
        ukvr = kvr1(np)
        uh = h1(np)
        uxmol(1:option%nspec) = xmol1((np-1)*option%nspec+1:np*option%nspec)
      else
        ukvr = kvr2(np)
        uh = h2(np)
        uxmol(1:option%nspec) = xmol2((np-1)*option%nspec+1:np*option%nspec)
      endif      
     
   ! print *,'FLcont  uxmol',uxmol
      if (ukvr>floweps) then
        v_darcy= Dq * ukvr * dphi
         !option%vvl_loc(nconn_no) = v_darcy
        vv_darcy(np) = v_darcy
     
        q = v_darcy * area
          
        do m=1, option%nspec 
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
      do m=1, option%nspec
        ind=m+(np-1)*option%nspec
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
  Res_FL(1:option%ndof-1) = fluxm(:) * option%dt
  Res_FL(option%ndof) = fluxe * option%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 !                           2 R = R - Res_FL  
 !print *,'end FLcont'
 
  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)       
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)       

end subroutine RichardsRes_FLCont

#ifdef OVERHAUL
#ifdef OVERHAUL_V2
subroutine RichardsRes_FLBCCont(nbc_no,ibc,area,var_node1,var_node2,por2,tor2,sir2, &
                              dd1,perm2,Dk2,dist_gravity,option,vv_darcy, &
                              Res_FL)
  use Option_module
  use Grid_module
#else
subroutine RichardsRes_FLBCCont(nbc_no,area,var_node1,var_node2,por2,tor2,sir2, &
                              dd1,perm2,Dk2,grid,vv_darcy,connection_object,Res_FL)
#endif                              
  use Connection_module

#else
subroutine RichardsRes_FLBCCont(nbc_no,area,var_node1,var_node2,por2,tor2,sir2, &
                              dd1,perm2,Dk2,grid,vv_darcy,Res_FL)
#endif
 ! Notice : index 1 stands for BC node
 
  implicit none
  
  integer nbc_no, ibc
  type(option_type) :: option
  real*8 dd1, sir2(1:option%nphase)
  real*8, target:: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8, target:: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  real*8 por2,perm2,Dk2,tor2
  real*8 vv_darcy(option%nphase), area
  real*8 Res_FL(1:option%ndof) 
  
  real*8 :: dist_gravity  ! distance along gravity vector
          
  real*8, pointer :: temp1, pre_ref1   ! 1 dof
  real*8, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  real*8, pointer :: xmol1(:), diff1(:)            ! 
  
  real*8, pointer :: temp2, pre_ref2   ! 1 dof
  real*8, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  real*8, pointer :: xmol2(:), diff2(:)    
  
  integer ibase, m,np, ind, j
  real*8  fluxm(option%nspec),fluxe, v_darcy,q
  real*8 uh,uxmol(1:option%nspec), ukvr,diff,diffdp, DK,Dq
  real*8 upweight,density_ave,cond,gravity, dphi

#ifndef OVERHAUL_V2
#ifdef OVERHAUL
    type(connection_type) :: connection_object
#endif
#endif
  
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
  ibase=ibase+option%nphase; 
                         xmol1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)
                         xmol2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)
  ibase=ibase+option%nphase*option%nspec;
                         diff1=>var_node1(ibase:ibase+option%nphase*option%nspec-1)    
                         diff2=>var_node2(ibase:ibase+option%nphase*option%nspec-1)


  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy = 0.D0 
   
  select case (option%ibndtyp(ibc))
    case(1) 
      Dq = perm2 / dd1
      diffdp = por2*tor2/dd1*area
        ! Flow term
      do np=1,option%nphase
        if ((sat1(np) > sir2(np)) .or. (sat2(np) > sir2(np))) then
          upweight=1.D0
        if (sat1(np) <eps) then 
          upweight=0.d0
        else if (sat2(np) <eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  

#ifndef OVERHAUL   
        gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * option%gravity * grid%delzbc(nbc_no)
#else
#ifdef OVERHAUL_V2
        gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * option%gravity * dist_gravity
#else
        gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * option%gravity * connection_object%dist(3,nbc_no) * &
                  connection_object%dist(0,nbc_no)
#endif       
#endif       
        dphi = pre_ref1- pre_ref2 + gravity
   
   
        if (dphi>=0.D0) then
          ukvr = kvr1(np)
          uh = h1(np)
          uxmol(:)=xmol1((np-1)*option%nspec+1 : np*option%nspec)
        else
          ukvr = kvr2(np)
          uh = h2(np)
          uxmol(:)=xmol2((np-1)*option%nspec+1 : np*option%nspec)
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        !option%vvl_loc(nbc_no) = v_darcy
          vv_darcy(np) = v_darcy
     
          q = v_darcy * area
          
          do m=1, option%nspec 
            fluxm(m) = fluxm(m) + q*density_ave*uxmol(m)
          enddo
          
          fluxe = fluxe + q*density_ave*uh 
        endif
      endif 
! Diffusion term   
! Note : average rule may not be correct  
      if ((sat1(np) > eps) .and. (sat2(np) > eps)) then
     
        diff = diffdp * 0.25D0*(sat1(np)+sat2(np))*(density1(np)+density2(np))
        do m = 1, option%nspec
          ind=m+(np-1)*option%nspec
          fluxm(m) = fluxm(m) + diff * diff2(ind)*( xmol1(ind)-xmol2(ind))
        enddo  
      endif
    enddo
! conduction term
        
    Dk =  Dk2 / dd1
    cond = Dk*area*(temp1-temp2) 
    fluxe=fluxe + cond
 
    Res_FL(1:option%nspec)=fluxm(:)* option%dt
    Res_FL(option%ndof)=fluxe * option%dt

  case(2)
    if ((dabs(option%velocitybc(1,nbc_no)))>floweps) then
!geh      print *, 'FlowBC :', nbc_no,option%velocitybc(1,nbc_no), &
!geh               option%velocitybc(2,nbc_no)
      
      do j=1,option%nphase
        v_darcy = option%velocitybc(j,nbc_no)
        vv_darcy(j) = option%velocitybc(j,nbc_no)
!      option%vvbc(j+(nc-2)*option%nphase)= option%velocitybc(j,nc)
      ! note different from 2 phase version

        if (v_darcy >0.d0) then 
          q = v_darcy * density1(j) * area
             !q = 0.d0
             !flux = flux - q
          fluxe = fluxe + q  * h1(j) 
          do m=1, option%nspec
            fluxm(m) = fluxm(m) + q * xmol1(m + (j-1)*option%nspec)
          enddo 
        else 
          q = v_darcy * density2(j) * area   
          fluxe = fluxe + q  * h2(j) 
          do m=1, option%nspec
            fluxm(m) = fluxm(m) + q * xmol2(m + (j-1)*option%nspec)
          enddo 
        endif 
   
      enddo
    endif
    Res_FL(1:option%nspec) = fluxm(:)*option%dt
    Res_FL(option%ndof) = fluxe * option%dt

  case(3)
    Dq = perm2/dd1 

    do np =1,option%nphase
      if ((sat2(np) > sir2(np))) then
    
        density_ave = density2(np)  
#ifndef OVERHAUL    
        gravity = density2(np)*amw2(np)* option%gravity * grid%delzbc(nbc_no)
#else
#ifdef OVERHAUL_V2
        gravity = density2(np)*amw2(np)* option%gravity * dist_gravity
#else
        gravity = density2(np)*amw2(np)* option%gravity * &
                  connection_object%dist(3,nbc_no) * &
                  connection_object%dist(0,nbc_no)
#endif
#endif        
        dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
    
        ukvr = kvr2(np)
        uh = h2(np)
        uxmol(:) = xmol2((np-1)*option%nspec+1 : np*option%nspec)
         
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
    !     option%vvbc(,nbc_no) = v_darcy
          vv_darcy(np) = v_darcy
     
          q = v_darcy * area
          
          do m=1, option%nspec 
            fluxm(m) = fluxm(m) + q*density_ave*uxmol(m)
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

end  subroutine RichardsRes_FLBCCont 


subroutine RichardsResidual(snes,xx,r,solution,ierr)

  use water_eos_module
  use Gas_Eos_Module
  use translator_Richards_module

#ifdef OVERHAUL  
  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
#endif  
  
  implicit none

#include "definitions.h"
 
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(solution_type) :: solution

! integer :: j, jm1, jm2, jmu, mu
  integer :: ierr
  integer :: n, ng, nc, nr
  integer :: i, i1, i2, jn, jng
  integer :: m, m1, m2, n1, n2, ip1, ip2, p1, p2
  !, t1, t2, c1, c2, s1, s2
  integer :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
! integer :: i1_hencoeff, i2_hencoeff
  integer :: ibc  ! Index that specifies a boundary condition block
  
! real*8 :: term1, term2, term3

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  PetscScalar, pointer ::accum_p(:)

  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
!              ddensity_p(:), ddensity_loc_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), var_p(:),var_loc_p(:) 
                          
               
! PetscScalar, pointer :: pc_p(:), pc_loc_p(:), kvr_p(:), kvr_loc_p(:)

  PetscScalar, pointer :: iphase_loc_p(:),icap_p(:),iphase_p(:),&
                          icap_loc_p(:), ithrm_loc_p(:),ithrm_p(:)

  integer :: iicap,iiphase, index_var_begin, index_var_end,iicap1,iicap2,np

  real*8 :: dd1, dd2, &
!           eng, cond, den, eengl,eengg, &
!           fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl,&
!           fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q,  &
!           v_darcy,hflx, &
            pvoldt, voldt, accum, pvol
  real*8 :: dd, f1, f2, ff
! real*8 :: por1, por2, density_ave
  real*8 :: perm1, perm2
! real*8 :: Dphi,D0
! real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
! real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol,dif(solution%option%nphase)
  real*8 :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1!, qqsrc
! real*8 :: cw,cw1,cw2, xxlw,xxla,xxgw,xxga
! real*8 :: upweight
! real*8 :: ukvr,uhh,uconc
  real*8 :: tmp
! real*8 :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp, visc
  real*8 :: rho
  real*8 :: Res(solution%option%ndof), vv_darcy(solution%option%nphase)
 PetscViewer :: viewer

#ifdef OVERHAUL 
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity
  
  grid => solution%grid
  option => solution%option
#endif
 
  option%vvlbc=0.D0
  option%vvgbc=0.D0
  option%vvl_loc=0.D0
  option%vvg_loc=0.D0

 
  call GridGlobalToLocal(grid,xx,option%xx_loc,NDOF)
  call GridGlobalToLocal(grid,option%iphas,option%iphas_loc,ONEDOF)

  call VecGetArrayF90(option%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%iphas_loc, iphase_loc_p, ierr); CHKERRQ(ierr)

! there is potential possiblity that the pertubation of p may change the direction of pflow.
! once that happens, code may crash, namely wrong derive. 
  do ng = 1, grid%ngmax 
    !ng=grid%nL2G(n)
    iiphase=int(iphase_loc_p(ng))
  
    if(xx_loc_p((ng-1)*option%ndof+1)>=option%pref)then
      option%delx(1,ng)=xx_loc_p((ng-1)*option%ndof+1)*dfac
    else
      option%delx(1,ng)=-xx_loc_p((ng-1)*option%ndof+1)*dfac
    endif  
      
    option%delx(2,ng)=xx_loc_p((ng-1)*option%ndof+2)*dfac
    if (option%ndof > 2) &
      option%delx(3:option%ndof,ng)=xx_loc_p((ng-1)*option%ndof+3:ng*option%ndof)*dfac
  
   enddo
  
  call VecRestoreArrayF90(option%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(option%iphas_loc, iphase_loc_p, ierr)

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%icap,icap_p,ierr)
  call VecGetArrayF90(option%ithrm,ithrm_p,ierr)  
  call VecGetArrayF90(option%iphas, iphase_p, ierr)
  call VecGetArrayF90(option%var,var_p,ierr)
  
  
 ! call VecGetArrayF90(option%ithrm,ithrm_p,ierr)
!------------------------------------------------------ 





!-----  phase properities ---- last time step---
  do n = 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    jn = 1 + (n-1)*option%nphase
    ng = grid%nL2G(n)
    ii1 = jn  !1+(n-1)*option%nphase
    ii2 = n*option%nphase
    iicap = icap_p(n)
    iiphase = iphase_p(n)
    !*****************
    dif(1)= option%difaq
    !dif(2)= option%cdiff(int(ithrm_p(n)))
  
  !*******************************************
    call pri_var_trans_richards_ninc(xx_p((n-1)*option%ndof+1:n*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use), &
                                option%itable,ierr, option%pref)

    iphase_p(n) =iiphase

    if (option%ideriv .eq. 1) then
      call pri_var_trans_richards_winc(xx_p((n-1)*option%ndof+1:n*option%ndof), &
                                  option%delx(1:option%ndof,ng),iiphase, &
                                  option%scale,option%nphase,option%nspec, &
                                  iicap ,dif, &
                                  var_p((n-1)*size_var_node+size_var_use+1:n* &
                                    size_var_node), &
                                  option%itable,ierr, option%pref)
    endif
  
!  print *,'var_p',n,iiphase, var_p((n-1)*size_var_node+1:n*size_var_node)              
!   if(n < 5) print *,'pflow_2ph: ',n,option%iideriv,option%xxphi_co2(n)
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(option%iphas,iphase_p,ierr)
  call VecRestoreArrayF90(option%icap,icap_p,ierr)
  call VecRestoreArrayF90(option%ithrm,ithrm_p,ierr)
  call VecRestoreArrayF90(option%var,var_p,ierr)
  ! call VecRestoreArrayF90(option%iphase,iphase_p,ierr)
  
  call GridGlobalToLocal(grid,option%var,option%var_loc,VARDOF)
  call GridGlobalToLocal(grid,option%perm_xx,option%perm_xx_loc,ONEDOF)
  call GridGlobalToLocal(grid,option%perm_yy,option%perm_yy_loc,ONEDOF)
  call GridGlobalToLocal(grid,option%perm_zz,option%perm_zz_loc,ONEDOF)
  call GridGlobalToLocal(grid,option%ithrm,option%ithrm_loc,ONEDOF)
  call GridGlobalToLocal(grid,option%icap,option%icap_loc,ONEDOF)

! End distribute data 
! now assign access pointer to local variables
  call VecGetArrayF90(option%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(option%accum, accum_p, ierr)
! call VecGetArrayF90(option%yy, yy_p, ierr)
 

  ! notice:: here we assume porosity is constant
 
  call VecGetArrayF90(option%var_loc,var_loc_p,ierr)
  call VecGetArrayF90(option%yy,yy_p,ierr)
  call VecGetArrayF90(option%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(option%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(option%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(option%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(option%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(option%volume, volume_p, ierr)
  call VecGetArrayF90(option%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(option%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(option%vl, vl_p, ierr)
  call VecGetArrayF90(option%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'


  if (option%rk > 0.d0) then
    call VecGetArrayF90(option%phis,phis_p,ierr)
  endif

  Resold_AR=0.D0; ResOld_FL=0.D0

!--------------------------------------------------------------------------
! Calculate accumulation term for interior and exterior nodes.
!--------------------------------------------------------------------------
! print *,option%rtot
  
!  print *, 'Residual  (init):'
!  print *, r_p
  
  r_p = - accum_p

!  print *, 'Residual  (after accum_p):'
!  print *, r_p

  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*option%ndof
    index_var_begin=(ng-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    
    pvol = volume_p(n)*porosity_loc_p(ng)
    voldt = volume_p(n) / option%dt
    pvoldt = porosity_loc_p(ng) * voldt
    iiphase = iphase_loc_p(ng)
    i = ithrm_loc_p(ng)

    accum = 0.d0
    call RichardsRes_ARCont(n, var_loc_p(index_var_begin: index_var_end),&
    porosity_loc_p(ng),volume_p(n),option%dencpr(i), option, Res, 1,ierr)
   
    r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
    Resold_AR(n,1:option%ndof)= Res(1:option%ndof) 
  enddo

!  print *, 'Residual  (after accum):'
!  print *, r_p

!************************************************************************
! add source/sink terms

!GEH - Structured Grid Dependence - Begin
  do nr = 1, option%nblksrc
      
    kk1 = option%k1src(nr) - grid%structured_grid%nzs
    kk2 = option%k2src(nr) - grid%structured_grid%nzs
    jj1 = option%j1src(nr) - grid%structured_grid%nys
    jj2 = option%j2src(nr) - grid%structured_grid%nys
    ii1 = option%i1src(nr) - grid%structured_grid%nxs
    ii2 = option%i2src(nr) - grid%structured_grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%structured_grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%structured_grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%structured_grid%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle

!geh - there has to be a better way than the below.
    do i = 2, option%ntimsrc
      if (option%timesrc(i,nr) == option%t) then
        tsrc1 = option%tempsrc(i,nr)
        qsrc1 = option%qsrc(i,nr)
        csrc1 = option%csrc(i,nr)
        hsrc1=  option%hsrc(i,nr)
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
    
   !print *,'pflow2ph : ', grid%myrank,i,option%timesrc(i,nr), &
   !option%timesrc(i-1,nr),option%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
   if(dabs(hsrc1)>1D-20)then 
       do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%structured_grid%nlx+(kk-1)*grid%structured_grid%nlxy
             r_p(n*option%ndof) = r_p(n*option%ndof) - hsrc1 * option%dt   
           enddo
          enddo
       enddo
  endif         

    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%structured_grid%nlx+(kk-1)*grid%structured_grid%nlxy
            ng = grid%nL2G(n)

            call wateos_noderiv(tsrc1,var_loc_p((ng-1)*size_var_node+2), &
                                dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
              
            r_p((n-1)*option%ndof + option%jh2o) = r_p((n-1)*option%ndof + option%jh2o) &
                                               - qsrc1 *option%dt
            r_p(n*option%ndof) = r_p(n*option%ndof) - qsrc1*enth_src_h2o*option%dt
            Resold_AR(n,option%jh2o)= Resold_AR(n,option%jh2o) - qsrc1*option%dt
            Resold_AR(n,option%ndof)= Resold_AR(n,option%ndof) - qsrc1 * &
                                                          enth_src_h2o * option%dt
      
      
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
            n = ii+(jj-1)*grid%structured_grid%nlx+(kk-1)*grid%structured_grid%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*option%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            call ideal_gaseos_noderiv(var_loc_p((ng-1)*size_var_node+2), &
                 tsrc1,option%scale,rho,enth_src_co2, tmp)
            enth_src_co2 =enth_src_co2 / option%fmwa

            r_p((n-1)*option%ndof + option%jco2) = r_p((n-1)*option%ndof + option%jco2) &
                                               - csrc1*option%dt
            r_p(n*option%ndof) = r_p(n*option%ndof) - csrc1 * enth_src_co2 *option%dt
            Resold_AR(n,option%jco2)= Resold_AR(n,option%jco2) - csrc1*option%dt
            Resold_AR(n,option%ndof)= Resold_AR(n,option%ndof) - csrc1 * &
                                                        enth_src_co2 * option%dt
       !r_p(s1) = r_p(s1) - csrc1

        !   print *,'pflow2ph_co2: ',grid%myrank,nr,n,ng,tsrc1,rho,option%fmwco2,csrc1
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

!GEH - Structured Grid Dependence - End

  !print *,'finished source/sink term'
  
!  print *, 'Residual  (after source/sink):'
!  print *, r_p

!*********************************************************************


 
! stop
!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------

#ifndef OVERHAUL
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)
#else
  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      m1 = cur_connection_object%id_up(iconn)
      m2 = cur_connection_object%id_dn(iconn)
      nc = iconn
#endif

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   

    if (associated(option%imat)) then
      if (option%imat(m1) <= 0 .or. option%imat(m2) <= 0) cycle
    endif

    p1 = 1 + (n1-1)*option%ndof 
    p2 = 1 + (n2-1)*option%ndof

#ifndef OVERHAUL    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End
#else
    fraction_upwind = cur_connection_object%dist(-1,iconn)
    distance = cur_connection_object%dist(0,iconn)
    ! The below assumes a unit gravity vector of [0,0,1]
    distance_gravity = cur_connection_object%dist(3,iconn)*distance
    dd1 = distance*fraction_upwind
    dd2 = distance-dd1 ! should avoid truncation error
    
    ! for now, just assume diagonal tensor
    perm1 = perm_xx_loc_p(m1)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(m1)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(m1)*cur_connection_object%dist(3,iconn)

    perm2 = perm_xx_loc_p(m2)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(m2)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(m2)*cur_connection_object%dist(3,iconn)
#endif

    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    iicap1=int(icap_loc_p(m1))
    iicap2=int(icap_loc_p(m2))
   
    D1 = option%ckwet(i1)
    D2 = option%ckwet(i2)

#ifndef OVERHAUL
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
#endif
!   if(dabs(perm1-1D-15)>1D-20)print *, 'perm1 error', perm1, ip1, n1,n2
!  if(dabs(perm2-1D-15)>1D-20)print *, 'perm2 error', perm2, ip2, n1,n2
#ifdef OVERHAUL
#ifdef OVERHAUL_V2
    call RichardsRes_FLCont(iconn,cur_connection_object%area(iconn), &
                          var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m1),tor_loc_p(m1), &
                          option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                          var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m2),tor_loc_p(m2), &
                          option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                          distance_gravity,fraction_upwind,option, &
                          vv_darcy,Res)
#else
    call RichardsRes_FLCont(iconn,cur_connection_object%area(iconn), &
                          var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m1),tor_loc_p(m1), &
                          option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                          var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m2),tor_loc_p(m2), &
                          option%sir(1:option%nphase,iicap2),dd2,perm2,D2,grid, &
                          vv_darcy,cur_connection_object,Res)
#endif
#else
    call RichardsRes_FLCont(nc ,option%area(nc), &
                          var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m1),tor_loc_p(m1), &
                          option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                          var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                            size_var_node+size_var_use), &
                          porosity_loc_p(m2),tor_loc_p(m2), &
                          option%sir(1:option%nphase,iicap2),dd2,perm2,D2,grid, &
                          vv_darcy,Res)
#endif
#ifndef OVERHAUL
    option%vvl_loc(nc) = vv_darcy(1)
!    option%vvg_loc(nc) = vv_darcy(2)  
#else
    cur_connection_object%velocity(1,iconn) = vv_darcy(1)
#endif    
    if (n1 > 0) then               ! If the upstream node is not a ghost node...
      do np =1, option%nphase 
#ifndef OVERHAUL      
        vl_p(np+(ip1-1)*option%nphase+3*option%nphase*(n1-1)) = vv_darcy(np) 
#else
        vl_p(np+(0)*option%nphase+3*option%nphase*(n1-1)) = &
                                     vv_darcy(np)*cur_connection_object%dist(1,iconn) 
        vl_p(np+(1)*option%nphase+3*option%nphase*(n1-1)) = &
                                     vv_darcy(np)*cur_connection_object%dist(2,iconn) 
        vl_p(np+(2)*option%nphase+3*option%nphase*(n1-1)) = &
                                     vv_darcy(np)*cur_connection_object%dist(3,iconn) 
#endif        
        ! use for print out of velocity
      enddo
    endif
     
    Resold_FL(nc,1:option%ndof) = Res(1:option%ndof) 
    
    if (n1>0) then
      r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
    endif
   
    if (n2>0) then
      r_p(p2:p2+option%ndof-1) = r_p(p2:p2+option%ndof-1) - Res(1:option%ndof)
    endif

#ifndef OVERHAUL
  enddo
#else
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
#endif
   ! print *,'finished NC' 
 
!  print *, 'Residual  (after flux):'
!  print *, r_p

!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', grid%myrank, option%sgbc    
 
#ifndef OVERHAUL
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
#else
  connection_list => grid%boundary_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      m= cur_connection_object%id_dn(iconn)
      nc = iconn
#endif
    ng = grid%nL2G(m)

    if (associated(option%imat)) then
      if (option%imat(ng) <= 0) cycle
    endif

    if (ng<=0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif

    p1 = 1 + (m-1) * option%ndof

    ibc = grid%ibconn(nc)
    
    i2 = ithrm_loc_p(ng)
    D2 = option%ckwet(i2)

#ifndef OVERHAUL    
!GEH - Structured Grid Dependence - Begin
    ip1 = grid%ipermbc(nc)
!GEH - Structured Grid Dependence - End

!GEH - Structured Grid Dependence - Begin
    select case(ip1)
      case(1)
        perm1 = perm_xx_loc_p(ng)
      case(2)
        perm1 = perm_yy_loc_p(ng)
      case(3)
        perm1 = perm_zz_loc_p(ng)
    end select
!GEH - Structured Grid Dependence - End

#else
    ! for now, just assume diagonal tensor
    perm1 = perm_xx_loc_p(ng)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(ng)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(ng)*cur_connection_object%dist(3,iconn)
    ! The below assumes a unit gravity vector of [0,0,1]
    distance_gravity = cur_connection_object%dist(3,iconn) * &
                       cur_connection_object%dist(0,iconn)
#endif

    select case(option%ibndtyp(ibc))
          
      case(2)
      ! solve for pb from Darcy's law given qb /= 0
!        option%xxbc(:,nc) = xx_loc_p((ng-1)*option%ndof+1: ng*option%ndof)
!        option%iphasebc(nc) = int(iphase_loc_p(ng))
        option%xxbc(:,nc)=xx_loc_p((ng-1)*option%ndof+1:ng*option%ndof)
         option%iphasebc(nc) = int(iphase_loc_p(ng))
        if(dabs(option%velocitybc(1,nc))>1D-20)then
          if( option%velocitybc(1,nc)>0) then
             option%xxbc(2:option%ndof,nc)= yybc(2:option%ndof,nc)
          endif     
        endif    


      case(3) 
     !  option%xxbc((nc-1)*option%ndof+1)=grid%pressurebc(2,ibc)
        option%xxbc(2:option%ndof,nc) = xx_loc_p((ng-1)*option%ndof+2: ng*option%ndof)
        option%iphasebc(nc)=int(iphase_loc_p(ng))
    
     case(4)
        option%xxbc(1,nc) = xx_loc_p((ng-1)*option%ndof+1)
         option%xxbc(3:option%ndof,nc) = xx_loc_p((ng-1)*option%ndof+3: ng*option%ndof)    
        option%iphasebc(nc)=int(iphase_loc_p(ng))
        
    end select

! print *,'2ph bc',grid%myrank,nc,m,ng,ibc,option%ibndtyp(ibc),grid%pressurebc(:,ibc), &
! option%tempbc(ibc),option%sgbc(ibc),option%concbc(ibc),option%velocitybc(:,ibc)

!   if(option%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   grid%ibndphaseRate(ibc) 
   !   grid%ibndconc(ibc)    ! 
   !   option%tempbc(ibc)      !1 elements 
   !   grid%pressurebc(ibc)  !nphase elements
!      endif
   
    
    iicap=int(icap_loc_p(ng))  
!      print *,'pflow_2pha_bc: ',grid%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      grid%nconnbc,option%ibndtyp(ibc),option%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,option%iideriv,option%ibndtyp(ibc),option%density_bc,&
!   grid%pressurebc(2,ibc),option%tempbc(ibc),option%concbc(ibc),option%sgbc(ibc)
       
        !*****************
    dif(1)= option%difaq
!    dif(2)= option%cdiff(int(ithrm_loc_p(ng)))
    !*******************************************

!    print *, nc
!    print *, 'xxbc: ', option%xxbc(:,nc)
!    print *, 'iphasebc: ', option%iphasebc(nc)
!    print *, 'icaptype: ', option%icaptype(iicap)
!    print *, 'sir: ', option%sir(1:option%nphase,iicap)
!    print *, 'lambda: ', grid%lambda(iicap)
!    print *, 'alpha: ', option%alpha(iicap)
!    print *, 'pckrm: ', grid%pckrm(iicap)
!    print *, 'pwrprm: ', grid%pwrprm(iicap)
!    print *, 'varbc: ', option%varbc(1:size_var_use)
!    print *, 'xxphi_co2_bc: ', option%xxphi_co2_bc(nc)
!    print *
  
    call pri_var_trans_Richards_ninc(option%xxbc(:,nc),option%iphasebc(nc),&
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif,&
                                option%varbc(1:size_var_use),option%itable,ierr, &
                                option%pref)
   
#ifdef OVERHAUL
#ifdef OVERHAUL_V2
    call RichardsRes_FLBCCont(nc,grid%ibconn(nc), &
                            cur_connection_object%area(iconn), &
                            option%varbc(1:size_var_use), &
                            var_loc_p((ng-1)*size_var_node+1:(ng-1)* &
                            size_var_node+size_var_use),porosity_loc_p(ng), &
                            tor_loc_p(ng),option%sir(1:option%nphase,iicap), &
                            cur_connection_object%dist(0,iconn),perm1,D2, &
                            distance_gravity,option, &
                            vv_darcy,Res)
#else
    call RichardsRes_FLBCCont(nc,cur_connection_object%area(iconn), &
                            option%varbc(1:size_var_use), &
                            var_loc_p((ng-1)*size_var_node+1:(ng-1)* &
                            size_var_node+size_var_use),porosity_loc_p(ng), &
                            tor_loc_p(ng),option%sir(1:option%nphase,iicap), &
                            cur_connection_object%dist(0,iconn),perm1,D2, grid, &
                            vv_darcy,cur_connection_object,Res)
#endif
#else
    call RichardsRes_FLBCCont(nc,option%areabc(nc),option%varbc(1:size_var_use), &
                            var_loc_p((ng-1)*size_var_node+1:(ng-1)* &
                            size_var_node+size_var_use),porosity_loc_p(ng), &
                            tor_loc_p(ng),option%sir(1:option%nphase,iicap), &
                            grid%distbc(nc),perm1,D2, grid, vv_darcy,Res)
#endif     
#ifndef OVERHAUL
    option%vvlbc(nc) = vv_darcy(1)
 !   option%vvgbc(nc) = vv_darcy(2) 
#else
    cur_connection_object%velocity(1,iconn) = vv_darcy(1)
#endif    
    r_p(p1:p1-1+option%ndof)= r_p(p1:p1-1+option%ndof) - Res(1:option%ndof)
    ResOld_AR(m,1:option%ndof) = ResOld_AR(m,1:option%ndof) - Res(1:option%ndof)
   
   
       !print *, ' boundary index', nc,ng,ibc,option%ibndtyp(ibc)
       !print *,'        xxbc', option%iphasebc(nc), option%xxbc(:,nc),res
     !print *, '       var', option%varbc
  !   print *, ' P  T   C   S  ', grid%pressurebc(1,ibc),option%tempbc(ibc), &
    !                               option%concbc(ibc),option%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),option%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,option%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),option%tempbc(ibc),option%concbc(ibc),option%sgbc(ibc)
!print *,option%density_bc,option%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc
#ifndef OVERHAUL
  enddo
#else
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
#endif
!  print *,'finished BC'

!  print *, 'Residual  (after bc flux):'
!  print *, r_p

  if (option%use_isoth==PETSC_TRUE) then
    do n = 1, grid%nlmax  ! For each local node do...

      !geh - Ignore inactive cells with inactive materials
      if (associated(option%imat)) then
        if (option%imat(grid%nL2G(n)) <= 0) cycle
      endif

      ng = grid%nL2G(n)   ! corresponding ghost index
      p1 = 3 + (n-1)*option%ndof
      r_p(p1)=xx_loc_p(2 + (ng-1)*option%ndof)-yy_p(p1-1)
    enddo
  endif

  if (n_zero_rows > 0) then
    do n=1,n_zero_rows
      r_p(zero_rows_local(n)) = 0.
    enddo
  endif

!  print *, 'Residual  (final):'
!  print *, r_p

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(option%yy, yy_p, ierr)
  call VecRestoreArrayF90(option%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(option%accum, accum_p, ierr)
  call VecRestoreArrayF90(option%var_loc,var_loc_p,ierr)
  call VecRestoreArrayF90(option%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(option%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(option%volume, volume_p, ierr)
  call VecRestoreArrayF90(option%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(option%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(option%vl, vl_p, ierr)
  call VecRestoreArrayF90(option%iphas_loc, iphase_loc_p, ierr)
  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(option%phis,phis_p,ierr)
  endif
#define DEBUG_GEH
#define DEBUG_GEH_ALL
#ifdef DEBUG_GEH 
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'residual.out',viewer,ierr)
 call VecView(r,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'xx.out',viewer,ierr)
 call VecView(xx,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

  !print *,'XX ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
 !print *,'Residual ::...........'; call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)
 
 !print *,'finished RichardsResidual'
end subroutine RichardsResidual
                
! --------------------------------------------------------------------- 

subroutine RichardsJacobian(snes,xx,A,B,flag,solution,ierr)
       
  use water_eos_module
  use gas_eos_module
  use translator_Richards_module

#ifdef OVERHAUL  
  use Connection_module
  use Option_module
  use Grid_module
  use Solution_module
#endif  
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(solution_type) :: solution
   ! integer, intent(inout) :: flag
  MatStructure flag

! integer :: j, jn, jm1, jm2, jmu, mu
  integer :: ierr
  integer :: n, ng, nc,nvar,neq,nr
  integer :: i1, i2, jng, i
  integer :: kk,ii1,jj1,kk1,ii2,jj2,kk2  
  integer :: m, m1, m2, n1, n2, ip1, ip2 
  integer :: p1,p2 !,t1,t2,c1,c2,s1,s2
  integer :: ibc  ! Index that specifies a boundary condition block.
! real*8 :: v_darcy, q
! real*8 :: dum1, dum2

  PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), phis_p(:),  tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  integer :: iicap,ii,jj,iiphas,iiphas1,iiphas2,iicap1,iicap2
  integer :: index_var_begin, index_var_end
! integer ibc_hencoeff
! real*8 :: dddt,dddp,fg,dfgdp
  real*8 :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
!           dfgdt,eng,dhdt,dhdp,visc,dvdt,dvdp
! real*8 :: cond, gravity, acc, density_ave, den
  real*8 :: vv_darcy(solution%option%nphase),voldt,pvoldt
! real*8 :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
!           flux, fluxh, fluxv, difff, diffg, diffl,
  real*8 :: ff,dif(1:solution%option%nphase)
  real*8 :: tsrc1,qsrc1,csrc1
  real*8 :: dd1, dd2, dd, f1, f2
! real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
! real*8 :: por1, por2
  real*8 :: perm1, perm2
! real*8 :: qu_rate, p_vapor,sat_pressure_t
! real*8 :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
! real*8 :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! real*8 :: sat_pressure  ! Saturation pressure of water.
! real*8 :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu, sat_ave
  real*8 :: ra(1:solution%option%ndof,1:2*solution%option%ndof)  
! real*8 :: uhh, uconc, ukvr
  real*8 :: tmp
! real*8 :: upweight,m1weight,m2weight,mbweight,mnweight
  real*8 :: delxbc(1:solution%option%ndof)
  real*8 :: blkmat11(1:solution%option%ndof,1:solution%option%ndof), &
            blkmat12(1:solution%option%ndof,1:solution%option%ndof),&
            blkmat21(1:solution%option%ndof,1:solution%option%ndof),&
            blkmat22(1:solution%option%ndof,1:solution%option%ndof)
  real*8 :: ResInc(1:solution%grid%nlmax, 1:solution%option%ndof, 1:solution%option%ndof),res(1:solution%option%ndof)  
  real*8 :: max_dev  
  integer  na1,na2
  
#ifdef OVERHAUL 
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option 
#endif  
  
 PetscViewer :: viewer
!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  grid => solution%grid
  option => solution%option

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
  flag = SAME_NONZERO_PATTERN

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
 ! call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
 !                           option%xx_loc, ierr)
 ! call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
 !                         option%xx_loc, ierr)
 
 !call DAGlobalToLocalBegin(grid%da_1_dof, option%iphas, &
 !                          INSERT_VALUES, option%iphas_loc, ierr)
 !call DAGlobalToLocalEnd(grid%da_1_dof, option%iphas, &
 !                         INSERT_VALUES, option%iphas_loc, ierr)

  call VecGetArrayF90(option%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(option%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(option%tor_loc, tor_loc_p, ierr)
! call VecGetArrayF90(option%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(option%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(option%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(option%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(option%volume, volume_p, ierr)

  call VecGetArrayF90(option%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(option%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(option%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(option%var_loc, var_loc_p, ierr)

 !print *,' In mph Jacobian ::  got pointers '
! ********************************************************************

! Accumulation terms

  ResInc=0.D0
  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ng = grid%nL2G(n)   !get ghosted index
    
    voldt = volume_p(n) / option%dt
    pvoldt = porosity_loc_p(ng) * voldt

    iiphas=iphase_loc_p(ng)
 ! pressure equation    
    do nvar=1, option%ndof
   
      index_var_begin=(ng-1)*size_var_node+nvar*size_var_use+1
      index_var_end = index_var_begin -1 + size_var_use

      call RichardsRes_ARCont(n, var_loc_p(index_var_begin : index_var_end), &
                            porosity_loc_p(ng),volume_p(n), &
                            option%dencpr(int(ithrm_loc_p(ng))),option, Res,1,ierr)
      
      ResInc(n,:,nvar) = ResInc(n,:,nvar) + Res(:)
    enddo
  enddo
! print *,' Mph Jaco Finished accum terms'
! Source / Sink term
#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

!GEH - Structured Grid Dependence - Begin
  do nr = 1, option%nblksrc
      
    kk1 = option%k1src(nr) - grid%structured_grid%nzs
    kk2 = option%k2src(nr) - grid%structured_grid%nzs
    jj1 = option%j1src(nr) - grid%structured_grid%nys
    jj2 = option%j2src(nr) - grid%structured_grid%nys
    ii1 = option%i1src(nr) - grid%structured_grid%nxs
    ii2 = option%i2src(nr) - grid%structured_grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%structured_grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%structured_grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%structured_grid%nlx,ii2)
        
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
10  continue
    
   !print *,'pflow2ph : ', grid%myrank,i,option%timesrc(i,nr), &
   !option%timesrc(i-1,nr),option%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
    
    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%structured_grid%nlx+(kk-1)*grid%structured_grid%nlxy
            ng = grid%nL2G(n)
            
            do nvar=1,option%ndof      
              call wateos_noderiv(tsrc1,var_loc_p((ng-1)*size_var_node+nvar* &
                                                         size_var_use+2), &
                                  dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
              ResInc(n,option%jh2o,nvar) = ResInc(n,option%jh2o,nvar) - qsrc1 * &
                                                                    option%dt
              ResInc(n,option%ndof,nvar) = ResInc(n,option%ndof,nvar) - qsrc1 * &
                                                        enth_src_h2o * option%dt

      
      
      !           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*option%fmwh2o, &
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
            n = ii+(jj-1)*grid%structured_grid%nlx+(kk-1)*grid%structured_grid%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*option%nphase
                    
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            do nvar=1,option%ndof     
              call ideal_gaseos_noderiv(var_loc_p((ng-1)*size_var_node+nvar* &
                                                         size_var_use+2), &
                                        tsrc1,option%scale,rho,enth_src_co2,tmp)
              enth_src_co2 = enth_src_co2 / option%fmwa
    
              ResInc(n,option%jco2,nvar)=  ResInc(n,option%jco2,nvar) - csrc1 * &
                                                                    option%dt
              ResInc(n,option%ndof,nvar)=  ResInc(n,option%ndof,nvar) - csrc1 * &
                                                           enth_src_co2*option%dt

          !  Res_AR(n,option%jco2)= Res_AR(n,option%jco2) - csrc1
    !  Res_AR(n,option%ndof)= Res_AR(n,option%ndof) - csrc1 * enth_src_co2
       !r_p(s1) = r_p(s1) - csrc1

!           print *,'pflow2ph_co2: ',nr,n,ng,tsrc1,rho,option%fmwco2,csrc1
            enddo
          enddo
        enddo
      enddo
    endif
  enddo  

!GEH - Structured Grid Dependence - End

  ! print *,' Mph Jaco Finished source terms'
! Contribution from BC
#ifndef OVERHAUL
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
#else
  connection_list => grid%boundary_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      m = cur_connection_object%id_dn(iconn)
      nc = iconn
#endif
    ng = grid%nL2G(m)

    if (associated(option%imat)) then
      if (option%imat(ng) <= 0) cycle
    endif

    if (ng<=0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif
  
    p1 = 1 + (m-1) * option%ndof
       
    ibc = grid%ibconn(nc)
    
    i2 = ithrm_loc_p(ng)
    D2 = option%ckwet(i2)

#ifndef OVERHAUL 
!GEH - Structured Grid Dependence - Begin
    ip1 = grid%ipermbc(nc)
!GEH - Structured Grid Dependence - End

!GEH - Structured Grid Dependence - Begin
    select case(ip1)
      case(1)
        perm1 = perm_xx_loc_p(ng)
      case(2)
        perm1 = perm_yy_loc_p(ng)
      case(3)
        perm1 = perm_zz_loc_p(ng)
    end select
!GEH - Structured Grid Dependence - End
#else
    ! for now, just assume diagonal tensor
    perm1 = perm_xx_loc_p(ng)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(ng)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(ng)*cur_connection_object%dist(3,iconn)
    ! The below assumes a unit gravity vector of [0,0,1]
    distance_gravity = cur_connection_object%dist(3,iconn) * &
                       cur_connection_object%dist(0,iconn)
#endif
    delxbc=0.D0
    select case(option%ibndtyp(ibc))
      case(1)
        delxbc =0.D0
      case(2)
        ! solve for pb from Darcy's law given qb /= 0
!        option%xxbc(:,nc) = xx_loc_p((ng-1)*option%ndof+1: ng*option%ndof)
!        option%iphasebc(nc) = int(iphase_loc_p(ng))
!        delxbc = option%delx(1:option%ndof,ng)
        option%xxbc(:,nc) = xx_loc_p((ng-1)*option%ndof+1: ng*option%ndof)
        option%iphasebc(nc) = int(iphase_loc_p(ng))
        delxbc = option%delx(1:option%ndof,ng)
        
  
       if(dabs(option%velocitybc(1,nc))>1D-20)then
         if( option%velocitybc(1,nc)>0) then
             option%xxbc(2:option%ndof,nc)= yybc(2:option%ndof,nc)
             delxbc(2:option%ndof)=0.D0
          endif     
        endif    

      
       
        
     case(3) 
        !    option%xxbc(1,nc)=grid%pressurebc(2,ibc)
        option%xxbc(2:option%ndof,nc) = xx_loc_p((ng-1)*option%ndof+2:ng*option%ndof)
        option%iphasebc(nc) = int(iphase_loc_p(ng))
        delxbc(1) = 0.D0
        delxbc(2:option%ndof) = option%delx(2:option%ndof,ng)
        
        case(4)
        option%xxbc(1,nc) = xx_loc_p((ng-1)*option%ndof+1)
        option%xxbc(3:option%ndof,nc) = xx_loc_p((ng-1)*option%ndof+3: ng*option%ndof)    
        delxbc(1)=option%delx(1,ng)
        delxbc(3:option%ndof) = option%delx(3:option%ndof,ng) 
        option%iphasebc(nc)=int(iphase_loc_p(ng))
    
    end select

! print *,'2ph bc',grid%myrank,nc,m,ng,ibc,option%ibndtyp(ibc),grid%pressurebc(:,ibc), &
! option%tempbc(ibc),option%sgbc(ibc),option%concbc(ibc),option%velocitybc(:,ibc)

!   if(option%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   grid%ibndphaseRate(ibc) 
   !   grid%ibndconc(ibc)    ! 
   !   option%tempbc(ibc)      !1 elements 
   !   grid%pressurebc(ibc)  !nphase elements
!      endif
   
    iicap = int(icap_loc_p(ng))     
       
!      print *,'pflow_2pha_bc: ',grid%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      grid%nconnbc,option%ibndtyp(ibc),option%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,option%iideriv,option%ibndtyp(ibc),option%density_bc,&
!   grid%pressurebc(2,ibc),option%tempbc(ibc),option%concbc(ibc),option%sgbc(ibc)
        !*****************
    dif(1) = option%difaq
  !  dif(2) = option%cdiff(int(ithrm_loc_p(ng)))
    !*******************************************

  !  print *,' Mph Jaco BC terms: finish setup'
  ! here should pay attention to BC type !!!
    call pri_var_trans_Richards_ninc(option%xxbc(:,nc),option%iphasebc(nc), &
                                option%scale,option%nphase,option%nspec, &
                                iicap,  dif, &
                                option%varbc(1:size_var_use),option%itable,ierr, option%pref)
  
    call pri_var_trans_Richards_winc(option%xxbc(:,nc),delxbc,option%iphasebc(nc), &
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif(1:option%nphase),&
                                option%varbc(size_var_use+1:(option%ndof+1)* &
                                  size_var_use), &
                                option%itable,ierr, option%pref)
            
!    print *,' Mph Jaco BC terms: finish increment'
    do nvar=1,option%ndof
#ifdef OVERHAUL
#ifdef OVERHAUL_V2
      call RichardsRes_FLBCCont(nc,grid%ibconn(nc), &
                              cur_connection_object%area(iconn), &
                              option%varbc(nvar*size_var_use+1:(nvar+1)* &
                                size_var_use), &
                              var_loc_p((ng-1)*size_var_node+nvar* &
                                size_var_use+1:(ng-1)*size_var_node+nvar* &
                                size_var_use+size_var_use), &
                              porosity_loc_p(ng),tor_loc_p(ng), &
                              option%sir(1:option%nphase,iicap), &
                              cur_connection_object%dist(0,iconn),perm1,D2, &
                              distance_gravity,option,vv_darcy, &
                              Res)
#else
      call RichardsRes_FLBCCont(nc,cur_connection_object%area(iconn), &
                              option%varbc(nvar*size_var_use+1:(nvar+1)* &
                                size_var_use), &
                              var_loc_p((ng-1)*size_var_node+nvar* &
                                size_var_use+1:(ng-1)*size_var_node+nvar* &
                                size_var_use+size_var_use), &
                              porosity_loc_p(ng),tor_loc_p(ng), &
                              option%sir(1:option%nphase,iicap), &
                              cur_connection_object%dist(0,iconn),perm1,D2, &
                              grid,vv_darcy,cur_connection_object,Res)
#endif                              
#else   
      call RichardsRes_FLBCCont(nc,option%areabc(nc), &
                              option%varbc(nvar*size_var_use+1:(nvar+1)* &
                                size_var_use), &
                              var_loc_p((ng-1)*size_var_node+nvar* &
                                size_var_use+1:(ng-1)*size_var_node+nvar* &
                                size_var_use+size_var_use), &
                              porosity_loc_p(ng),tor_loc_p(ng), &
                              option%sir(1:option%nphase,iicap), &
                              grid%distbc(nc),perm1,D2,grid,vv_darcy,Res)
#endif    
      ResInc(m,1:option%ndof,nvar) = ResInc(m,1:option%ndof,nvar) - Res(1:option%ndof)
    enddo
 !  print *,' Mph Jaco BC terms: finish comp'
    !   print *, ' boundary index', nc,ng,ibc,option%ibndtyp(ibc)
    !   print *, ' P  T   C   S  ', grid%pressurebc(1,ibc),option%tempbc(ibc), &
    !                               option%concbc(ibc),option%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),option%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,option%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),option%tempbc(ibc),option%concbc(ibc),option%sgbc(ibc)
!print *,option%density_bc,option%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc

#ifndef OVERHAUL
  enddo
#else
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
#endif
  ! print *,' Mph Jaco Finished BC terms'
#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
  do n= 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    ra=0.D0
    ng = grid%nL2G(n)
    na1= grid%nG2N(ng)
   ! Remember, the matrix index starts from (0,0)
    p1 = (ng-1)*option%ndof ! = 1 + (ng-1)*option%ndof-1
   
    max_dev = 0.D0
    do neq=1, option%ndof
      do nvar=1, option%ndof
        ra(neq,nvar) = ResInc(n,neq,nvar)/option%delx(nvar,ng) - &
                       ResOld_AR(n,neq)/option%delx(nvar,ng)
        if (max_dev < dabs(ra(neq,nvar))) max_dev = dabs(ra(neq,nvar))
   
      enddo      
    enddo
    if (option%use_isoth==PETSC_TRUE) then
      ra(:,2)=0.D0
      ra(3,1:option%ndof)=0.D0
      ra(3,2)=1.D0
    endif     
   
    if (max_dev<1D-5) then
      print *,'Mph Jaco max dev = ', max_dev
    endif
  
    if (option%iblkfmt == 0) then
      p1=(na1)*option%ndof
      do ii=0,option%ndof-1
        do jj=0,option%ndof-1
          call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
        enddo
      enddo
    else
      blkmat11=ra(1:option%ndof,1:option%ndof)
      call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
    endif
         
  enddo
#ifdef DEBUG_GEH_ALL    
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
  
!   print *,' Mph Jaco Finished one node terms'
! -----------------------------contribution from transport----------------------

 !print *,'phase cond: ',iphase_loc_p
  ResInc=0.D0
  
#ifndef OVERHAUL  
  do nc = 1, grid%nconn  ! For each interior connection...
    ra = 0.D0
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)
#else
  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      m1 = cur_connection_object%id_up(iconn)
      m2 = cur_connection_object%id_dn(iconn)
      nc = iconn
#endif

    if (associated(option%imat)) then
      if (option%imat(m1) <= 0 .or. option%imat(m2) <= 0) cycle
    endif

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   
    na1 = grid%nG2N(m1)
    na2 = grid%nG2N(m2)
    !print *, grid%myrank,nc,m1,m2,n1,n2,na1,na2
    p1 =  (m1-1)*option%ndof
    p2 =  (m2-1)*option%ndof
   
#ifndef OVERHAUL    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)

!GEH - Structured Grid Dependence - Begin    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)
!GEH - Structured Grid Dependence - End

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

#else
    fraction_upwind = cur_connection_object%dist(-1,iconn)
    distance = cur_connection_object%dist(0,iconn)
    ! The below assumes a unit gravity vector of [0,0,1]
    distance_gravity = cur_connection_object%dist(3,iconn)*distance
    dd1 = distance*fraction_upwind
    dd2 = distance-dd1 ! should avoid truncation error
    
    ! for now, just assume diagonal tensor
    perm1 = perm_xx_loc_p(m1)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(m1)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(m1)*cur_connection_object%dist(3,iconn)

    perm2 = perm_xx_loc_p(m2)*cur_connection_object%dist(1,iconn)+ &
            perm_yy_loc_p(m2)*cur_connection_object%dist(2,iconn)+ &
            perm_zz_loc_p(m2)*cur_connection_object%dist(3,iconn)
#endif
    
    iiphas1 = iphase_loc_p(m1)
    iiphas2 = iphase_loc_p(m2)

    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = option%ckwet(i1)
    D2 = option%ckwet(i2)

#ifndef OVERHAUL
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
#endif    
    
    iicap1 = int(icap_loc_p(m1))
    iicap2 = int(icap_loc_p(m2))
 
  ! do neq = 1, option%ndof
    do nvar = 1, option%ndof
#ifdef OVERHAUL    
#ifdef OVERHAUL_V2  
      call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                            var_loc_p((m1-1)*size_var_node+nvar* &
                              size_var_use+1:(m1-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                            distance_gravity,fraction_upwind, &
                            option,vv_darcy,Res)
#else
      call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                            var_loc_p((m1-1)*size_var_node+nvar* &
                              size_var_use+1:(m1-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                            grid,vv_darcy,cur_connection_object,Res)
#endif
#else
      call RichardsRes_FLCont(nc ,option%area(nc), &
                            var_loc_p((m1-1)*size_var_node+nvar* &
                              size_var_use+1:(m1-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+1:(m2-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2,grid, &
                            vv_darcy,Res)
#endif
      ra(:,nvar)= (Res(:)-ResOld_FL(nc,:))/option%delx(nvar,m1)
       
#ifdef OVERHAUL 
#ifdef OVERHAUL_V2   
      call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                            var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+nvar* &
                              size_var_use+1:(m2-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                            distance_gravity,fraction_upwind, &
                            option, &
                            vv_darcy,Res)
#else
      call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                            var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+nvar* &
                              size_var_use+1:(m2-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2,grid, &
                            vv_darcy,cur_connection_object,Res)
#endif
#else
      call RichardsRes_FLCont(nc,option%area(nc), &
                            var_loc_p((m1-1)*size_var_node+1:(m1-1)* &
                              size_var_node+size_var_use),&
                            porosity_loc_p(m1),tor_loc_p(m1), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((m2-1)*size_var_node+nvar* &
                              size_var_use+1:(m2-1)*size_var_node+nvar* &
                              size_var_use+size_var_use),&
                            porosity_loc_p(m2),tor_loc_p(m2), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2,grid, &
                            vv_darcy,Res)
#endif
 
      ra(:,nvar+option%ndof)= (Res(:)-ResOld_FL(nc,:))/option%delx(nvar,m2)
   
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
    p1=(na1)*option%ndof;p2=(na2)*option%ndof
    do ii=0,option%ndof-1
      do jj=0,option%ndof-1
        if (n1>0) then
          if (option%iblkfmt == 0) then
            call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
          else
            blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
          endif
        endif
        if (n2>0) then
          if (option%iblkfmt == 0) then
            call MatSetValue(A,p2+ii,p1+jj,-ra(ii+1,jj+1),ADD_VALUES,ierr)
          else
            blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
          endif
        endif
      enddo
   
      do jj=option%ndof,2*option%ndof-1
        if (n1>0) then
          if (option%iblkfmt == 0) then
            call MatSetValue(A,p1+ii,p2+jj-option%ndof,ra(ii+1,jj+1), &
                             ADD_VALUES,ierr)
          else
            blkmat12(ii+1,jj-option%ndof+1) = blkmat12(ii+1,jj-option%ndof+1) + &
                                                             ra(ii+1,jj+1)
          endif
        endif
        if (n2>0) then
          if (option%iblkfmt == 0) then
            call MatSetValue(A,p2+ii,p2+jj-option%ndof,-ra(ii+1,jj+1), &
                             ADD_VALUES,ierr)
          else
            blkmat22(ii+1,jj-option%ndof+1) =  blkmat22(ii+1,jj-option%ndof+1) - &
                                                              ra(ii+1,jj+1)
          endif
        endif
      enddo
    enddo
  
    if (option%iblkfmt /= 0) then
      if (n1>0) call MatSetValuesBlocked(A,1,na1,1,na1,blkmat11,ADD_VALUES,ierr)
      if (n2>0) call MatSetValuesBlocked(A,1,na2,1,na2,blkmat22,ADD_VALUES,ierr)
      if (n1>0) call MatSetValuesBlocked(A,1,na1,1,na2,blkmat12,ADD_VALUES,ierr)
      if (n2>0) call MatSetValuesBlocked(A,1,na2,1,na1,blkmat21,ADD_VALUES,ierr)
    endif
!print *,'accum r',ra(1:5,1:8)   
 !print *,'devq:',nc,q,dphi,devq(3,:)
#ifndef OVERHAUL
  enddo
#else
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
#endif
  ! print *,' Mph Jaco Finished Two node terms'
  
  call VecRestoreArrayF90(option%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(option%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(option%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(option%var_loc, var_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(option%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(option%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(option%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(option%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(option%iphas_loc, iphase_loc_p, ierr)

  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(option%phis,phis_p,ierr)
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
    if (n == 1) then
      p2 = p1
    elseif (n == 0) then
      p2 = p1-option%ndof+2
    else
      p2 = p1+1
    endif
    call MatSetValuesLocal(A,1,p1,1,p2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  f1 = 1.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,f1,ierr) 
#endif

  !B = A
  !call MatCopy(A,B,ierr)
  
 !call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB, ierr)

#ifdef DEBUG_GEH    
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine RichardsJacobian




subroutine pflow_Richards_initaccum(solution)
 
  use translator_Richards_module
  use Solution_module
  use Grid_module  
  use Option_module
 
  implicit none
  
  type(solution_type) :: solution 

 
  integer :: ierr
  integer :: n
  integer :: i, index_var_begin,index_var_end
  integer :: p1
! integer :: ii1,ii2
  integer :: iicap, iiphase


  PetscScalar, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_p(:),&
                          var_p(:), icap_p(:),iphase_p(:),ithrm_p(:)
  
 !  integer, pointer ::iphase_p(:)
  
! real*8 :: sat_pressure, pvol, satw  ! Saturation pressure of water.
  real*8 :: dif(1:solution%option%nphase),res(1:solution%option%ndof)
 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  grid => solution%grid
  option => solution%option
 
  call VecGetArrayF90(option%volume, volume_p, ierr)
  call VecGetArrayF90(option%porosity, porosity_p, ierr)
  call VecGetArrayF90(option%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%accum, accum_p, ierr)
  call VecGetArrayF90(option%var, var_p,ierr)
  call VecGetArrayF90(option%iphas, iphase_p, ierr)
  call VecGetArrayF90(option%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(option%icap, icap_p, ierr)
  !print *,'Richardsinitaccum  Gotten pointers'
 
  do n = 1, grid%nlmax
        
    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    iicap=int(icap_p(n))
    iiphase = int(iphase_p(n))
    dif(1)= option%difaq
 !   dif(2)= option%cdiff(int(ithrm_p(n)))

    call pri_var_trans_Richards_ninc(yy_p((n-1)*option%ndof+1:n*option%ndof),iiphase,&
                                option%scale,option%nphase,option%nspec, &
                                iicap , dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                size_var_node+size_var_use),option%itable,ierr, &
                                option%pref)
  enddo

  call VecRestoreArrayF90(option%var, var_p,ierr)
  call VecGetArrayF90(option%var, var_p,ierr)

!---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

  !  ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*option%ndof
    index_var_begin=(n-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    i = ithrm_p(n)
    
    call RichardsRes_ARCont(n,var_p(index_var_begin:index_var_end), &
                          porosity_p(n),volume_p(n),option%dencpr(i),option,Res, &
                          0,ierr)
 

    accum_p(p1:p1+option%ndof-1)=Res(:) 

   !print *, 'init m accum ', n,  Res 

! print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*option%nspec+(n-1)*option%nphase*option%nspec),kvr_p(jn),kvr_p(jn+1)

  enddo

  call VecRestoreArrayF90(option%volume, volume_p, ierr)
  call VecRestoreArrayF90(option%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(option%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(option%accum, accum_p, ierr)
  call VecRestoreArrayF90(option%var, var_p,ierr)
  call VecRestoreArrayF90(option%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(option%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(option%icap, icap_p, ierr)

end subroutine pflow_Richards_initaccum


subroutine pflow_update_Richards(solution)

  use translator_Richards_module
  use pckr_module
  use Condition_module_old
   ! use water_eos_module
#ifdef OVERHAUL  
  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
#endif 

  implicit none

  type(solution_type) :: solution 
    
! integer :: ichange 
  integer :: n,n0
!geh added for transient boundary conditons
  integer :: nc, ibc, iithrm, m, ng
  real*8 :: sw, pc(2), kr(2)
  integer :: ierr,iicap,iiphase, iiphase_old
  PetscScalar, pointer :: xx_p(:),icap_p(:),ithrm_p(:),iphase_p(:), var_p(:), yy_p(:), iphase_old_p(:)
  real*8 :: dif(1:solution%option%nphase)
! real*8 :: dum1, dum2           

#ifdef OVERHAUL 
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  grid => solution%grid
  option => solution%option
#endif
      
!geh added for transient boundary conditions      
  if (associated(option%imat) .and. option%iread_geom < 0) then
!commend out for now    call UpdateBoundaryConditions(option)
    yybc =option%xxbc
    vel_bc = option%velocitybc
  endif
!geh end

  ! if (option%rk > 0.d0) call Rock_Change(grid)
  ! call  Translator_Richards_Switching(option%xx,grid,1,ichange)
  !print *,'Richards_Update done'
 
   ! if(ichange ==1)then
  call VecGetArrayF90(option%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(option%icap,icap_p,ierr)
  call VecGetArrayF90(option%ithrm,ithrm_p,ierr)  
  call VecGetArrayF90(option%iphas, iphase_p, ierr)
  call VecGetArrayF90(option%iphas_old, iphase_old_p, ierr)
  call VecGetArrayF90(option%var,var_p,ierr)

  do n = 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    iicap = icap_p(n)
    iiphase = int(iphase_p(n))
    iiphase_old = int(iphase_old_p(n))
    n0 = (n-1)*option%ndof
    
     if(option%ndof>=3) then
        if (xx_p(n0+3)<0.D0) xx_p(n0+3)=1.D-6
     endif
    !*****************
  !  if (iiphase ==3)then
  !    if(xx_p(n0+1) > yy_p(n0+1) .and. xx_p(n0+1)+ formeps > option%pref .and. iiphase_old==3)then
  !     iphase_p(n) =1
  !     iiphase =1
  !     xx_p(n0+1) = option%pref + formeps
  !    endif
  !   elseif (iiphase ==1)then
  !    if(xx_p(n0+1) < yy_p(n0+1) .and. xx_p(n0+1) - formeps < option%pref .and. iiphase_old==1)then
  !     iphase_p(n) =3
  !     iiphase =3
  !     xx_p(n0+1) = option%pref - formeps
  !    endif
  !  endif   

    
    dif(1) = option%difaq
 !   dif(2) = option%cdiff(int(ithrm_p(n)))
    !*******************************************
    call pri_var_trans_Richards_ninc(xx_p((n-1)*option%ndof+1:n*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif,&
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use),&
                                option%itable,ierr, option%pref)
  ! print *,n, xx_p((n-1)*option%ndof+1:n*option%ndof), var_p((n-1)*size_var_node+1:(n-1)*
  !size_var_node+4)
 
   enddo

  !geh added for transient boundary conditions  
  if (associated(option%imat) .and. option%iread_geom < 0) then
#ifndef OVERHAUL
    do nc = 1, grid%nconnbc

      m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
#else
    connection_list => grid%boundary_connection_list
    cur_connection_object => connection_list%first
    do 
      if (.not.associated(cur_connection_object)) exit
      do iconn = 1, cur_connection_object%num_connections
        m = cur_connection_object%id_dn(iconn)
        nc = iconn
#endif

      ng = grid%nL2G(m)

      if (associated(option%imat)) then
        if (option%imat(ng) <= 0) cycle
      endif
       
      if (m<0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ibc = grid%ibconn(nc)

       
    !print *,'initadj_bc',nc,ibc,option%ibndtyp(ibc),grid%nconnbc

      if (option%ibndtyp(ibc)==1 .or.option%ibndtyp(ibc)==3) then
        iicap=int(icap_p(m))
        iithrm=int(ithrm_p(m)) 
        dif(1)= option%difaq
      !  dif(2)= option%cdiff(iithrm)
      
        if(option%iphasebc(nc) ==3)then
          sw= option%xxbc(1,nc)
          call pflow_pckr_richards_fw(iicap ,sw,pc,kr)    
          !if(pc(1)>grid%pcwmax(iicap))then
          !  print *,'INIT Warning: Pc>pcmax'
          !  pc(1)=grid%pcwmax(iicap)
          !endif 
          option%xxbc(1,nc) =  option%pref - pc(1)
        endif
      
        call pri_var_trans_Richards_ninc(option%xxbc(:,nc),option%iphasebc(nc), &
                                         option%scale,option%nphase,option%nspec, &
                                         iicap,dif, &
                                         option%varbc(1:size_var_use),option%itable,ierr, &
                                         option%pref)
      
        if (translator_check_cond_Richards(option%iphasebc(nc), &
                                           option%varbc(1:size_var_use), &
                                           option%nphase,option%nspec) /=1) then
          print *," Wrong bounday node init...  STOP!!!", option%xxbc(:,nc)
      
          print *,option%varbc
          stop    
        endif 
      endif

      if (option%ibndtyp(ibc)==2) then
        yybc(2:option%ndof,nc)= option%xxbc(2:option%ndof,nc)
        vel_bc(1,nc) = option%velocitybc(1,nc)
!        print *,'initadj', nc, yybc(:,nc), vel_bc(:,nc)
      endif 
      
#ifndef OVERHAUL
    enddo
#else
      enddo
      cur_connection_object => cur_connection_object%next
    enddo
#endif
  endif
 
  call VecRestoreArrayF90(option%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(option%yy, yy_p, ierr);
  call VecRestoreArrayF90(option%icap,icap_p,ierr)
  call VecRestoreArrayF90(option%ithrm,ithrm_p,ierr)  
  call VecRestoreArrayF90(option%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(option%iphas_old, iphase_old_p, ierr)
  call VecRestoreArrayF90(option%var,var_p,ierr)
   
  call translator_Richards_massbal(solution)
 ! endif 

  call VecCopy(option%xx, option%yy, ierr)   
  call VecCopy(option%iphas, option%iphas_old, ierr)   
   
  call  pflow_Richards_initaccum(solution)
    !print *,'pflow_Richards_initaccum done'
  call translator_Richards_get_output(grid%nlmax,option)
 ! print *,'translator_get_output done'
  ! the output variables should be put into grid%pressure, temp,xmol,sat...
  ! otherwise need to rewrite the pflow_output

end subroutine pflow_update_Richards





subroutine pflow_Richards_initadj(solution)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  use translator_Richards_module  
  use pckr_module, only: pflow_pckr_richards_fw
#ifdef OVERHAUL  
  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
#endif 
  
  implicit none

  type(solution_type) :: solution 

 
  integer :: ierr
  integer :: n, nc, ng
  integer :: ibc,jn
  integer :: m
  integer :: ii1,ii2,iicap
  integer :: iiphase,iithrm

  PetscScalar, pointer :: xx_p(:),var_p(:)
  PetscScalar, pointer ::iphase_p(:), ithrm_p(:),icap_p(:)
  
  real*8 :: dif(solution%option%nphase)
! real*8 :: dum1, dum2
  real*8 :: pc(1:solution%option%nphase), kr(1:solution%option%nphase), sw

#ifdef OVERHAUL 
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  grid => solution%grid
  option => solution%option  
#endif
  
! real*8 :: temp1
!  real*8, parameter :: Rg=8.31415D0

  call VecGetArrayF90(option%icap, icap_p, ierr)
  call VecGetArrayF90(option%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(option%iphas, iphase_p, ierr)
  call VecGetArrayF90(option%xx, xx_p, ierr)
  call VecGetArrayF90(option%var, var_p, ierr)
! print *,'initadj gotten pointers' 


  do n = 1, grid%nlmax
   
    !geh - Ignore inactive cells with inactive materials
    if (associated(option%imat)) then
      if (option%imat(grid%nL2G(n)) <= 0) cycle
    endif

    jn = 1 + (n-1)*option%nphase
    ii1=1+(n-1)*option%nphase; ii2=n*option%nphase
    iicap=int(icap_p(n))
        
    iiphase = iphase_p(n)
        !*****************
    dif(1)= option%difaq
 !   dif(2)= option%cdiff(int(ithrm_p(n)))
    !*******************************************
   if(iiphase ==3)then
    ! print *, 'rich iniadj: ',n, iiphase, xx_p((n-1)*option%ndof+1: n*option%ndof)
     sw= xx_p((n-1)*option%ndof+1)
     call pflow_pckr_richards_fw(iicap,sw,pc,kr)    
  !   print *,'INIT ', sw, pc(1), iicap, grid%pcwmax(iicap), option%icaptype(iicap),option%sir(1,iicap), grid%lambda(iicap), &
  !               option%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),sw,pc,kr,&
  !               grid%pcbetac(iicap),grid%pwrprm(iicap)
                 
     if(pc(1)>option%pcwmax(iicap))then
        print *,'INIT Warning: Pc>pcmax', sw, pc(1), iicap, option%pcwmax(iicap)
        pc(1)=option%pcwmax(iicap)
     endif 
     xx_p((n-1)*option%ndof+1)= option%pref - pc(1)
     print *,'Richards: Conv: ',n, sw, iicap, sw, pc(1),xx_p((n-1)*option%ndof+1: n*option%ndof)
    endif
    
    call pri_var_trans_Richards_ninc(xx_p((n-1)*option%ndof+1:n*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap,  dif, &
                                var_p((n-1)*size_var_node+1:(n-1)* &
                                  size_var_node+size_var_use), &
                                option%itable,ierr, option%pref)
   
    if (translator_check_cond_Richards(iiphase, &
                                        var_p((n-1)*size_var_node+1:(n-1)* &
                                        size_var_node+size_var_use), &
                                        option%nphase,option%nspec) /= 1 ) then
      print *," Wrong internal node init...  STOP!!!"
      stop    
    endif 
  enddo

  allocate(yybc(option%ndof,grid%boundary_connection_list%first%num_connections))
  allocate(vel_bc(option%nphase,grid%boundary_connection_list%first%num_connections))
  yybc =option%xxbc
  vel_bc = option%velocitybc

#ifndef OVERHAUL
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
#else
  connection_list => grid%boundary_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      m = cur_connection_object%id_dn(iconn)
      nc = iconn
#endif

    ng = grid%nL2G(m)

    if (associated(option%imat)) then
      if (option%imat(ng) <= 0) cycle
    endif
       
    if (m<0) then
      print *, "Wrong boundary node index... STOP!!!"
      stop
    endif

    ibc = grid%ibconn(nc)
       
!geh      print *,'initadj_bc',nc,ibc,option%ibndtyp(ibc),grid%nconnbc

    if (option%ibndtyp(ibc)==1 .or.option%ibndtyp(ibc)==3) then
      iicap=int(icap_p(m))
      iithrm=int(ithrm_p(m)) 
      dif(1)= option%difaq
    !  dif(2)= option%cdiff(iithrm)
      
     if(option%iphasebc(nc) ==3)then
       sw= option%xxbc(1,nc)
       call pflow_pckr_richards_fw(iicap,sw,pc,kr)    
    !   if(pc(1)>grid%pcwmax(iicap))then
    !     print *,'INIT Warning: Pc>pcmax'
    !     pc(1)=grid%pcwmax(iicap)
    !   endif 
        option%xxbc(1,nc) =  option%pref - pc(1)
     endif

      
      
      call pri_var_trans_Richards_ninc(option%xxbc(:,nc),option%iphasebc(nc), &
                                  option%scale,option%nphase,option%nspec, &
                                  iicap, dif, &
                                  option%varbc(1:size_var_use),option%itable,ierr, &
                                  option%pref)
      
      if (translator_check_cond_Richards(option%iphasebc(nc), &
                                          option%varbc(1:size_var_use), &
                                          option%nphase,option%nspec) /=1) then
        print *," Wrong bounday node init...  STOP!!!", option%xxbc(:,nc)
      
        print *,option%varbc
        stop    
      endif 
    endif

   if (option%ibndtyp(ibc)==2) then
  
       yybc(2:option%ndof,nc)= option%xxbc(2:option%ndof,nc)
       vel_bc(1,nc) = option%velocitybc(1,nc)
!geh       print *,'initadj', nc, yybc(:,nc), vel_bc(:,nc)
    endif 
#ifndef OVERHAUL
  enddo
#else
    enddo
    cur_connection_object => cur_connection_object%next
  enddo
#endif

  call VecRestoreArrayF90(option%icap, icap_p, ierr)
  call VecRestoreArrayF90(option%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(option%iphas, iphase_p, ierr)
  call VecRestoreArrayF90(option%xx, xx_p, ierr)
  call VecRestoreArrayF90(option%var, var_p, ierr)
  !print *,kgjkdf
  
  !call VecCopy(option%iphas,option%iphas_old,ierr)
   
end subroutine pflow_Richards_initadj


subroutine createRichardsZeroArray(solution)

  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none

  type(solution_type) :: solution
  integer :: n, ng, ncount, idof

#ifdef OVERHAUL
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  grid => solution%grid
  option => solution%option
#endif
  
  n_zero_rows = 0

  if (associated(option%imat)) then
    do n = 1, grid%nlmax
      ng = grid%nL2G(n)
      if (option%imat(ng) <= 0) then
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

  if (associated(option%imat)) then
    do n = 1, grid%nlmax
      ng = grid%nL2G(n)
      if (option%imat(ng) <= 0) then
        do idof = 1, option%ndof
          ncount = ncount + 1
          zero_rows_local(ncount) = (n-1)*option%ndof+idof
          zero_rows_local_ghosted(ncount) = (ng-1)*option%ndof+idof-1
        enddo
      else
#ifdef ISOTHERMAL
        ncount = ncount + 1
        zero_rows_local(ncount) = n*option%ndof
        zero_rows_local_ghosted(ncount) = ng*option%ndof-1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    do n = 1, grid%nlmax
      ng = grid%nL2G(n)
      ncount = ncount + 1
      zero_rows_local(ncount) = n*option%ndof
      zero_rows_local_ghosted(ncount) = ng*option%ndof-1
    enddo
#endif
  endif

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine createRichardsZeroArray


end module Richards_module
