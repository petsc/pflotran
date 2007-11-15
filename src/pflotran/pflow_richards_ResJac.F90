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
  use Field_module
 
  implicit none
  
  type(solution_type) :: solution
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  PetscScalar, pointer :: xx_p(:),yy_p(:)
  integer :: dof_offset,re,ierr
  integer :: local_id

  grid => solution%grid
  option => solution%option
  field => solution%field
 
  call VecGetArrayF90(field%xx, xx_p, ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr)

  do local_id=1, grid%nlmax
    dof_offset=(local_id-1)*option%ndof
    do re = 1, option%ndof
      xx_p(dof_offset+re)= yy_p(dof_offset+re)
    enddo
  enddo 
  call VecRestoreArrayF90(field%xx, xx_p, ierr) 
  call VecRestoreArrayF90(field%yy, yy_p, ierr)
 
end subroutine pflow_Richards_timecut
  

subroutine pflow_Richards_setupini(solution)

  use Solution_module
  use Option_module
  use Grid_module
  use Field_module
  use Region_module
  use Structured_Grid_module
  use Coupler_module
  use Condition_module
 
  implicit none
  
  type(solution_type) :: solution
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: initial_condition

  PetscScalar, pointer :: xx_p(:), iphase_loc_p(:)
  integer local_id, ghosted_id, ibegin, iend, icell, ierr
  
  grid => solution%grid
  option => solution%option
  field => solution%field
  
  size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
  size_var_node = (option%ndof + 1) * size_var_use
  
  allocate(Resold_AR(grid%nlmax,option%ndof))
  allocate(Resold_FL(grid%internal_connection_list%first%num_connections, &
                     option%ndof))
  allocate(option%delx(option%ndof,grid%ngmax))
  option%delx=0.D0
   
  call VecGetArrayF90(field%xx,xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  initial_condition => solution%initial_conditions%first
  
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
  
  call VecRestoreArrayF90(field%xx,xx_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)

end  subroutine pflow_Richards_setupini
  

subroutine Richards_Update_Reason(reason,solution)

  use Solution_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none
 
  type(solution_type) :: solution
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  integer, intent(out):: reason
  PetscScalar, pointer :: xx_p(:),iphase_loc_p(:), yy_p(:) !,r_p(:)
  integer :: dof_offset, temp_reason
  integer ierr, iipha
  integer :: local_id, ghosted_id
  
  grid => solution%grid
  option => solution%option
  field => solution%field
  
! real*8, pointer :: sat(:),xmol(:)
! real*8 rmax(option%ndof)

! why this barrier
  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  reason = 1
 ! call SNESComputeFunction(solver%snes,field%xx,grid%r,ierr)
 ! do n=1,option%ndof
 !  call VecStrideNorm(grid%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if (reason>0) then
    call VecGetArrayF90(field%xx,xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(field%yy,yy_p, ierr)
    call VecGetArrayF90(field%iphas_loc,iphase_loc_p, ierr); 
  
    do local_id = 1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      dof_offset=(local_id-1)* option%ndof 
      iipha=int(iphase_loc_p(ghosted_id))

      select case(iipha)
        case (1)
          if (xx_p(dof_offset + 1) < option%pref) then
            reason = 0
            exit
          endif
        case (3)
          if (xx_p(dof_offset + 1) > option%pref) then
            reason = 0
            exit
          endif
      end select  
    enddo
  
    if (reason<=0) print *,'Sat or Con out of Region at: ',local_id,iipha,xx_p(dof_offset+1:dof_offset+2)
    call VecRestoreArrayF90(field%xx,xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%yy,yy_p, ierr)
    call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p, ierr) 
  endif
 ! print *,' update reason', grid%myrank, re,n,grid%nlmax

  
  if (option%commsize >1) then
    temp_reason = reason
    call MPI_ALLREDUCE(temp_reason,reason,1, MPI_INTEGER,MPI_SUM, &
                       PETSC_COMM_WORLD,ierr)
    if (reason<option%commsize) reason = 0
  endif
  
  if (reason<=0) print *,'Sat or Con out of Region'
  
end subroutine Richards_Update_Reason

!=======================================================================================
 


subroutine RichardsRes_ARCont(node_no,var_node,por,vol,rock_dencpr,option,Res_AR,ireac,ierr)

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
  if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!H2O
    mol(1)= mol(1) - option%dt * option%rtot(node_no,1)
  endif
  Res_AR(1:option%ndof-1)=mol(:)
  Res_AR(option%ndof)=eng
  nullify(temp, pre_ref, sat, density, amw, h,u, pc,kvr,xmol,diff)       

end subroutine  RichardsRes_ARCont

subroutine RichardsRes_FLCont(nconn_no,area,var_node1,por1,tor1,sir1,dd1,perm1, &
                            Dk1,var_node2,por2,tor2,sir2,dd2,perm2,Dk2, &
                            dist_gravity,upweight, &
                            option,vv_darcy,Res_FL)
  use Option_module                              
  
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
   
  Dq = (perm1 * perm2)/(dd1*perm2 + dd2*perm1)
  diffdp = (por1 *tor1 * por2*tor2) / (dd2*por1*tor1 + dd1*por2*tor2)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy = 0.D0  
  
  do np=1, option%nphase

! Flow term
    if ((sat1(np) > sir1(np)) .or. (sat2(np) > sir2(np))) then
      if (sat1(np) <eps) then 
        upweight=0.d0
      else if (sat2(np) <eps) then 
        upweight=1.d0
      endif
      density_ave = upweight*density1(np)+(1.D0-upweight)*density2(np)  

      gravity = (upweight*density1(np)*amw1(np) + &
                (1.D0-upweight)*density2(np)*amw2(np)) &
                * option%gravity * dist_gravity

      dphi = pre_ref1 - pre_ref2  + gravity

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
     

      if (ukvr>floweps) then
        v_darcy= Dq * ukvr * dphi
        vv_darcy(np) = v_darcy
     
        q = v_darcy * area
          
        do m=1, option%nspec 
          fluxm(m)=fluxm(m) + q*density_ave*uxmol(m)
        enddo

        fluxe = fluxe + q*density_ave*uh 
      endif
    endif 

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

  enddo
   
! conduction term
        
  Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)
  cond = Dk*area*(temp1-temp2) 
  fluxe=fluxe + cond

  Res_FL(1:option%ndof-1) = fluxm(:) * option%dt
  Res_FL(option%ndof) = fluxe * option%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 !                           2 R = R - Res_FL  

 
  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)       
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)       

end subroutine RichardsRes_FLCont

subroutine RichardsRes_FLBCCont(nbc_no,ibndtype,area,var_node1,var_node2,por2,tor2,sir2, &
                              dd1,perm2,Dk2,dist_gravity,option,field,vv_darcy, &
                              Res_FL)
  use Option_module
  use Field_module
 
  implicit none
  
  integer nbc_no, ibndtype
  type(option_type) :: option
  type(field_type) :: field
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
   
  select case (ibndtype)
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

        gravity = (upweight*density1(np)*amw1(np) + &
                  (1.D0-upweight)*density2(np)*amw2(np)) &
                  * option%gravity * dist_gravity
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
    if ((dabs(field%velocitybc(1,nbc_no)))>floweps) then
      
      do j=1,option%nphase
        v_darcy = field%velocitybc(j,nbc_no)
        vv_darcy(j) = field%velocitybc(j,nbc_no)
      ! note different from 2 phase version

        if (v_darcy >0.d0) then 
          q = v_darcy * density1(j) * area
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
        gravity = density2(np)*amw2(np)* option%gravity * dist_gravity
        dphi = pre_ref1-pc1(np) - pre_ref2 + pc2(np) + gravity
    
        ukvr = kvr2(np)
        uh = h2(np)
        uxmol(:) = xmol2((np-1)*option%nspec+1 : np*option%nspec)
         
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
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

  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
  use Coupler_module  
  
  implicit none

#include "definitions.h"
 
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(solution_type) :: solution

  integer :: ierr
  integer :: nc
  integer :: i, ithrm1, ithrm2, jn
  integer :: ip1, ip2, p1, p2
  integer :: ii1,ii2
  integer :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  

  PetscScalar, pointer ::accum_p(:)

  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), var_p(:),var_loc_p(:) 
                          
               
  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  integer :: iicap,iiphase, index_var_begin, index_var_end,iicap1,iicap2,np

  real*8 :: dd1, dd2, &
            pvoldt, voldt, accum, pvol
  real*8 :: dd, f1, f2, ff
  real*8 :: perm1, perm2
  real*8 :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
  real*8 :: dw_kg, dw_mol,dif(solution%option%nphase)
  real*8 :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  real*8 :: tmp, upweight
  real*8 :: rho
  real*8 :: Res(solution%option%ndof), vv_darcy(solution%option%nphase)
 PetscViewer :: viewer

  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  logical :: enthalpy_flag
  integer :: iconn
  integer :: sum_connection
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity
  
  grid => solution%grid
  option => solution%option
  field => solution%field

  call GridGlobalToLocal(grid,xx,field%xx_loc,NDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)

  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); CHKERRQ(ierr)

! there is potential possiblity that the pertubation of p may change the direction of pflow.
! once that happens, code may crash, namely wrong derive. 
  do ghosted_id = 1, grid%ngmax 
    iiphase=int(iphase_loc_p(ghosted_id))
  
    if(xx_loc_p((ghosted_id-1)*option%ndof+1)>=option%pref)then
      option%delx(1,ghosted_id)=xx_loc_p((ghosted_id-1)*option%ndof+1)*dfac
    else
      option%delx(1,ghosted_id)=-xx_loc_p((ghosted_id-1)*option%ndof+1)*dfac
    endif  
      
    option%delx(2,ghosted_id)=xx_loc_p((ghosted_id-1)*option%ndof+2)*dfac
    if (option%ndof > 2) &
      option%delx(3:option%ndof,ghosted_id)=xx_loc_p((ghosted_id-1)*option%ndof+3:ghosted_id*option%ndof)*dfac
  
   enddo
  
  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%var_loc,var_loc_p,ierr)
  
!-----  phase properities ---- last time step---
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    jn = 1 + (local_id-1)*option%nphase
    ii1 = jn 
    ii2 = local_id*option%nphase
    iicap = icap_loc_p(ghosted_id)
    iiphase = iphase_loc_p(ghosted_id)

    dif(1)= option%difaq
  
  !*******************************************
    call pri_var_trans_richards_ninc(xx_p((local_id-1)*option%ndof+1: &
                                       local_id*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif, &
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                  size_var_node+size_var_use), &
                                option%itable,ierr, option%pref)

    iphase_loc_p(ghosted_id) = iiphase

    if (option%ideriv .eq. 1) then
      call pri_var_trans_richards_winc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof), &
                                  option%delx(1:option%ndof,ghosted_id),iiphase, &
                                  option%scale,option%nphase,option%nspec, &
                                  iicap ,dif, &
                                  var_loc_p((ghosted_id-1)*size_var_node+size_var_use+1:ghosted_id* &
                                    size_var_node), &
                                  option%itable,ierr, option%pref)
    endif

  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecRestoreArrayF90(field%var_loc,var_loc_p,ierr)
  
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
  
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    p1 = 1 + (local_id-1)*option%ndof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin-1 + size_var_use
    
    pvol = volume_p(local_id)*porosity_loc_p(ghosted_id)
    voldt = volume_p(local_id) / option%dt
    pvoldt = porosity_loc_p(ghosted_id) * voldt
    iiphase = iphase_loc_p(ghosted_id)
    i = ithrm_loc_p(ghosted_id)

    accum = 0.d0
    call RichardsRes_ARCont(local_id,var_loc_p(index_var_begin:index_var_end),&
                            porosity_loc_p(ghosted_id),volume_p(local_id), &
                            option%dencpr(i),option,Res,1,ierr)
   
    r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
    Resold_AR(local_id,1:option%ndof)= Res(1:option%ndof) 
  enddo

!************************************************************************
! add source/sink terms
  source_sink => solution%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (size(source_sink%condition%cur_value) > RICHARDS_CONCENTRATION_DOF) then
      enthalpy_flag = .true.
    else
      enthalpy_flag = .false.
    endif

    qsrc1 = source_sink%condition%cur_value(RICHARDS_PRESSURE_DOF)
    tsrc1 = source_sink%condition%cur_value(RICHARDS_TEMPERATURE_DOF)
    csrc1 = source_sink%condition%cur_value(RICHARDS_CONCENTRATION_DOF)
    if (enthalpy_flag) hsrc1 = source_sink%condition%cur_value(RICHARDS_ENTHALPY_DOF)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_object => source_sink%connection
    
    do iconn = 1, cur_connection_object%num_connections      
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (enthalpy_flag) then
        r_p(local_id*option%ndof) = r_p(local_id*option%ndof) - hsrc1 * option%dt   
      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+2), &
                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
              
        r_p((local_id-1)*option%ndof + option%jh2o) = r_p((local_id-1)*option%ndof + option%jh2o) &
                                               - qsrc1 *option%dt
        r_p(local_id*option%ndof) = r_p(local_id*option%ndof) - qsrc1*enth_src_h2o*option%dt
        Resold_AR(local_id,option%jh2o)= Resold_AR(local_id,option%jh2o) - qsrc1*option%dt
        Resold_AR(local_id,option%ndof)= Resold_AR(local_id,option%ndof) - qsrc1 * &
                                                             enth_src_h2o * option%dt
      endif  
    
      if (csrc1 > 0.d0) then ! injection
#if 0   
! this if for co2, which is not supported in richards   
        jng= 2 + (ng-1)*option%nphase    
 
     !  span-wagner
        call ideal_gaseos_noderiv(var_loc_p((ghosted_id-1)*size_var_node+2), &
                                  tsrc1,option%scale,rho,enth_src_co2, tmp)
        enth_src_co2 =enth_src_co2 / option%fmwa

            r_p((local_id-1)*option%ndof + option%jco2) = r_p((local_id-1)*option%ndof + option%jco2) &
                                               - csrc1*option%dt
            r_p(local_id*option%ndof) = r_p(local_id*option%ndof) - csrc1 * enth_src_co2 *option%dt
            Resold_AR(local_id,option%jco2)= Resold_AR(local_id,option%jco2) - csrc1*option%dt
            Resold_AR(local_id,option%ndof)= Resold_AR(local_id,option%ndof) - csrc1 * &
                                                        enth_src_co2 * option%dt
       !r_p(s1) = r_p(s1) - csrc1

        !   print *,'pflow2ph_co2: ',grid%myrank,nr,n,ng,tsrc1,rho,option%fmwco2,csrc1
#endif
      endif
  
  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
      
  !  endif
    enddo
    source_sink => source_sink%next
  enddo

!*********************************************************************

!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------

  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      nc = sum_connection

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

      call RichardsRes_FLCont(iconn,cur_connection_object%area(iconn), &
                              var_loc_p((ghosted_id_up-1)*size_var_node+1:(ghosted_id_up-1)* &
                                size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                              option%sir(1:option%nphase,iicap1),dd1,perm1,D1, &
                              var_loc_p((ghosted_id_dn-1)*size_var_node+1:(ghosted_id_dn-1)* &
                                size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                              option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                              distance_gravity,upweight,option, &
                              vv_darcy,Res)

      cur_connection_object%velocity(1,iconn) = vv_darcy(1)

      if (local_id_up > 0) then               ! If the upstream node is not a ghost node...
        do np =1, option%nphase 
          vl_p(np+(0)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                                       vv_darcy(np)*cur_connection_object%dist(1,iconn) 
          vl_p(np+(1)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                                       vv_darcy(np)*cur_connection_object%dist(2,iconn) 
          vl_p(np+(2)*option%nphase+3*option%nphase*(local_id_up-1)) = &
                                       vv_darcy(np)*cur_connection_object%dist(3,iconn) 
          ! use for print out of velocity
        enddo
      endif
     
      Resold_FL(nc,1:option%ndof) = Res(1:option%ndof) 
    
      if (local_id_up>0) then
        r_p(p1:p1+option%ndof-1) = r_p(p1:p1+option%ndof-1) + Res(1:option%ndof)
      endif
   
      if (local_id_dn>0) then
        r_p(p2:p2+option%ndof-1) = r_p(p2:p2+option%ndof-1) - Res(1:option%ndof)
      endif

    enddo
    cur_connection_object => cur_connection_object%next
  enddo    
 
  boundary_condition => solution%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection
    
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      nc = sum_connection
    
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      p1 = 1 + (local_id-1) * option%ndof

    
      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_object%dist(3,iconn))
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = abs(cur_connection_object%dist(3,iconn)) * &
                         cur_connection_object%dist(0,iconn)

      select case(boundary_condition%condition%itype(1))
          
        case(2)
        ! solve for pb from Darcy's law given qb /= 0
          field%xxbc(:,nc)=xx_loc_p((ghosted_id-1)*option%ndof+1:ghosted_id*option%ndof)
           field%iphasebc(nc) = int(iphase_loc_p(ghosted_id))
          if(dabs(field%velocitybc(1,nc))>1D-20)then
            if( field%velocitybc(1,nc)>0) then
               field%xxbc(2:option%ndof,nc)= yybc(2:option%ndof,nc)
            endif     
          endif    


        case(3) 
          field%xxbc(2:option%ndof,nc) = xx_loc_p((ghosted_id-1)*option%ndof+2:ghosted_id*option%ndof)
          field%iphasebc(nc)=int(iphase_loc_p(ghosted_id))
    
       case(4)
          field%xxbc(1,nc) = xx_loc_p((ghosted_id-1)*option%ndof+1)
           field%xxbc(3:option%ndof,nc) = xx_loc_p((ghosted_id-1)*option%ndof+3:ghosted_id*option%ndof)    
          field%iphasebc(nc)=int(iphase_loc_p(ghosted_id))
        
      end select 
    
      iicap=int(icap_loc_p(ghosted_id))  

      dif(1)= option%difaq
  
      call pri_var_trans_Richards_ninc(field%xxbc(:,nc),field%iphasebc(nc),&
                                       option%scale,option%nphase,option%nspec, &
                                       iicap, dif,&
                                       field%varbc(1:size_var_use),option%itable,ierr, &
                                       option%pref)
   
      call RichardsRes_FLBCCont(nc,boundary_condition%condition%itype(1), &
                                cur_connection_object%area(iconn), &
                                field%varbc(1:size_var_use), &
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                  size_var_node+size_var_use),porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id),option%sir(1:option%nphase,iicap), &
                                cur_connection_object%dist(0,iconn),perm1,D2, &
                                distance_gravity,option,field, &
                                vv_darcy,Res)
      cur_connection_object%velocity(1,iconn) = vv_darcy(1)

      r_p(p1:p1-1+option%ndof)= r_p(p1:p1-1+option%ndof) - Res(1:option%ndof)
      ResOld_AR(local_id,1:option%ndof) = ResOld_AR(local_id,1:option%ndof) - Res(1:option%ndof)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  if (option%use_isoth==PETSC_TRUE) then
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

end subroutine RichardsResidual
                
! --------------------------------------------------------------------- 

subroutine RichardsJacobian(snes,xx,A,B,flag,solution,ierr)
       
  use water_eos_module
  use gas_eos_module
  use translator_Richards_module

  use Connection_module
  use Option_module
  use Grid_module
  use Solution_module
  use Coupler_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(solution_type) :: solution
  MatStructure flag

  integer :: ierr
  integer :: nc,nvar,neq,nr
  integer :: ithrm1, ithrm2, i
  integer :: ip1, ip2 
  integer :: p1,p2

  PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), phis_p(:),  tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  integer :: iicap,iiphas,iiphas1,iiphas2,iicap1,iicap2
  integer :: ii, jj
  integer :: index_var_begin, index_var_end
  real*8 :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  real*8 :: vv_darcy(solution%option%nphase),voldt,pvoldt
  real*8 :: ff,dif(1:solution%option%nphase)
  real*8 :: tsrc1,qsrc1,csrc1,hsrc1
  real*8 :: dd1, dd2, dd, f1, f2
  real*8 :: perm1, perm2
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.

  real*8 :: ra(1:solution%option%ndof,1:2*solution%option%ndof)  
  real*8 :: tmp, upweight
  real*8 :: delxbc(1:solution%option%ndof)
  real*8 :: blkmat11(1:solution%option%ndof,1:solution%option%ndof), &
            blkmat12(1:solution%option%ndof,1:solution%option%ndof),&
            blkmat21(1:solution%option%ndof,1:solution%option%ndof),&
            blkmat22(1:solution%option%ndof,1:solution%option%ndof)
  real*8 :: ResInc(1:solution%grid%nlmax, 1:solution%option%ndof, 1:solution%option%ndof),res(1:solution%option%ndof)  
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
  real*8 :: distance_gravity 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option 
  type(field_type), pointer :: field  
  
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
  field => solution%field

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
  flag = SAME_NONZERO_PATTERN

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p, ierr)

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

      call RichardsRes_ARCont(local_id,var_loc_p(index_var_begin:index_var_end), &
                            porosity_loc_p(ghosted_id),volume_p(local_id), &
                            option%dencpr(int(ithrm_loc_p(ghosted_id))),option, Res,1,ierr)
      
      ResInc(local_id,:,nvar) = ResInc(local_id,:,nvar) + Res(:)
    enddo
  enddo

! Source / Sink term
#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian1.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif

  source_sink => solution%source_sinks%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (size(source_sink%condition%cur_value) > RICHARDS_CONCENTRATION_DOF) then
      enthalpy_flag = .true.
    else
      enthalpy_flag = .false.
    endif

    qsrc1 = source_sink%condition%cur_value(RICHARDS_PRESSURE_DOF)
    tsrc1 = source_sink%condition%cur_value(RICHARDS_TEMPERATURE_DOF)
    csrc1 = source_sink%condition%cur_value(RICHARDS_CONCENTRATION_DOF)
    if (enthalpy_flag) hsrc1 = source_sink%condition%cur_value(RICHARDS_ENTHALPY_DOF)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_object => source_sink%connection
    
    do iconn = 1, cur_connection_object%num_connections      
      local_id= cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
          
      if (qsrc1 > 0.d0) then ! injection
     
        do nvar=1,option%ndof      
          call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+nvar* &
                                                     size_var_use+2), &
                              dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!       units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!       qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
          ResInc(local_id,option%jh2o,nvar) = ResInc(local_id,option%jh2o,nvar) - qsrc1 * &
                                                                    option%dt
          ResInc(local_id,option%ndof,nvar) = ResInc(local_id,option%ndof,nvar) - qsrc1 * &
                                                        enth_src_h2o * option%dt
        enddo
      endif  
    
      if (csrc1 > 0.d0) then ! injection
#if 0      
        jng= 2 + (ghosted_id-1)*option%nphase

         !  span-wagner
        do nvar=1,option%ndof     
          call ideal_gaseos_noderiv(var_loc_p((ghosted_id-1)*size_var_node+nvar* &
                                                     size_var_use+2), &
                                    tsrc1,option%scale,rho,enth_src_co2,tmp)
          enth_src_co2 = enth_src_co2 / option%fmwa
    
          ResInc(n,option%jco2,nvar)=  ResInc(n,option%jco2,nvar) - csrc1 * &
                                                                    option%dt
          ResInc(n,option%ndof,nvar)=  ResInc(n,option%ndof,nvar) - csrc1 * &
                                                           enth_src_co2*option%dt
        enddo
#endif        
      endif
    enddo
    source_sink => source_sink%next
  enddo  

! Contribution from BC
  boundary_condition => solution%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection

    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      nc = sum_connection
    
      local_id = cur_connection_object%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif
  
      p1 = 1 + (local_id-1) * option%ndof
       
    
      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_object%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_object%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_object%dist(3,iconn))
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = abs(cur_connection_object%dist(3,iconn)) * &
                         cur_connection_object%dist(0,iconn)

      delxbc=0.D0
      select case(boundary_condition%condition%itype(1))
        case(1)
          delxbc =0.D0
        case(2)
          ! solve for pb from Darcy's law given qb /= 0
          field%xxbc(:,nc) = xx_loc_p((ghosted_id-1)*option%ndof+1: ghosted_id*option%ndof)
          field%iphasebc(nc) = int(iphase_loc_p(ghosted_id))
          delxbc = option%delx(1:option%ndof,ghosted_id)
        
  
          if(dabs(field%velocitybc(1,nc))>1D-20)then
            if( field%velocitybc(1,nc)>0) then
             field%xxbc(2:option%ndof,nc)= yybc(2:option%ndof,nc)
             delxbc(2:option%ndof)=0.D0
            endif     
          endif    

        case(3) 
          field%xxbc(2:option%ndof,nc) = xx_loc_p((ghosted_id-1)*option%ndof+2:ghosted_id*option%ndof)
          field%iphasebc(nc) = int(iphase_loc_p(ghosted_id))
          delxbc(1) = 0.D0
          delxbc(2:option%ndof) = option%delx(2:option%ndof,ghosted_id)
        
        case(4)
          field%xxbc(1,nc) = xx_loc_p((ghosted_id-1)*option%ndof+1)
          field%xxbc(3:option%ndof,nc) = xx_loc_p((ghosted_id-1)*option%ndof+3: ghosted_id*option%ndof)    
          delxbc(1)=option%delx(1,ghosted_id)
          delxbc(3:option%ndof) = option%delx(3:option%ndof,ghosted_id) 
          field%iphasebc(nc)=int(iphase_loc_p(ghosted_id))
    
      end select

      iicap = int(icap_loc_p(ghosted_id))     
       
      dif(1) = option%difaq

  ! here should pay attention to BC type !!!
      call pri_var_trans_Richards_ninc(field%xxbc(:,nc),field%iphasebc(nc), &
                                  option%scale,option%nphase,option%nspec, &
                                  iicap,  dif, &
                                  field%varbc(1:size_var_use),option%itable,ierr, option%pref)
  
      call pri_var_trans_Richards_winc(field%xxbc(:,nc),delxbc,field%iphasebc(nc), &
                                  option%scale,option%nphase,option%nspec, &
                                  iicap, dif(1:option%nphase),&
                                  field%varbc(size_var_use+1:(option%ndof+1)* &
                                    size_var_use), &
                                  option%itable,ierr, option%pref)
              
      do nvar=1,option%ndof
        call RichardsRes_FLBCCont(nc,boundary_condition%condition%itype(1), &
                                cur_connection_object%area(iconn), &
                                field%varbc(nvar*size_var_use+1:(nvar+1)* &
                                  size_var_use), &
                                var_loc_p((ghosted_id-1)*size_var_node+nvar* &
                                  size_var_use+1:(ghosted_id-1)*size_var_node+nvar* &
                                  size_var_use+size_var_use), &
                                porosity_loc_p(ghosted_id),tor_loc_p(ghosted_id), &
                                option%sir(1:option%nphase,iicap), &
                                cur_connection_object%dist(0,iconn),perm1,D2, &
                                distance_gravity,option,field,vv_darcy, &
                                Res)

        ResInc(local_id,1:option%ndof,nvar) = ResInc(local_id,1:option%ndof,nvar) - Res(1:option%ndof)
      enddo

    enddo
    boundary_condition => boundary_condition%next
  enddo

#ifdef DEBUG_GEH_ALL  
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian2.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
  do local_id= 1, grid%nlmax

    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(grid%nL2G(local_id)) <= 0) cycle
    endif

    ra=0.D0
    ghosted_id = grid%nL2G(local_id)
    natural_id_up= grid%nG2N(ghosted_id)
   ! Remember, the matrix index starts from (0,0)
    p1 = (ghosted_id-1)*option%ndof 
   
    max_dev = 0.D0
    do neq=1, option%ndof
      do nvar=1, option%ndof
        ra(neq,nvar) = ResInc(local_id,neq,nvar)/option%delx(nvar,ghosted_id) - &
                       ResOld_AR(local_id,neq)/option%delx(nvar,ghosted_id)
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
      p1=(natural_id_up)*option%ndof
      do ii=0,option%ndof-1
        do jj=0,option%ndof-1
          call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
        enddo
      enddo
    else
      blkmat11=ra(1:option%ndof,1:option%ndof)
      call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_up,blkmat11,ADD_VALUES,ierr)
    endif
         
  enddo
#ifdef DEBUG_GEH_ALL    
 call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
 call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian3.out',viewer,ierr)
 call MatView(A,viewer,ierr)
 call PetscViewerDestroy(viewer,ierr)
#endif
  
! -----------------------------contribution from transport----------------------

  ResInc=0.D0
  
  connection_list => grid%internal_connection_list
  cur_connection_object => connection_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      nc = sum_connection
    
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
 
      do nvar = 1, option%ndof
        call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                              var_loc_p((ghosted_id_up-1)*size_var_node+nvar* &
                                size_var_use+1:(ghosted_id_up-1)*size_var_node+nvar* &
                                size_var_use+size_var_use),&
                              porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                              option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                              var_loc_p((ghosted_id_dn-1)*size_var_node+1:(ghosted_id_dn-1)* &
                                size_var_node+size_var_use),&
                              porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                              option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                              distance_gravity,upweight, &
                              option,vv_darcy,Res)
                              
        ra(:,nvar)= (Res(:)-ResOld_FL(nc,:))/option%delx(nvar,ghosted_id_up)
       
        call RichardsRes_FLCont(nc,cur_connection_object%area(iconn), &
                              var_loc_p((ghosted_id_up-1)*size_var_node+1:(ghosted_id_up-1)* &
                                size_var_node+size_var_use),&
                              porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                              option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                              var_loc_p((ghosted_id_dn-1)*size_var_node+nvar* &
                                size_var_use+1:(ghosted_id_dn-1)*size_var_node+nvar* &
                                size_var_use+size_var_use),&
                              porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                              option%sir(1:option%nphase,iicap2),dd2,perm2,D2, &
                              distance_gravity,upweight, &
                              option, &
                              vv_darcy,Res)
 
        ra(:,nvar+option%ndof)= (Res(:)-ResOld_FL(nc,:))/option%delx(nvar,ghosted_id_dn)
   
      enddo
   
      if (option%use_isoth==PETSC_TRUE) then
        ra(3,1:2*option%ndof)=0.D0
        ra(:,2)=0.D0
        ra(:,2+option%ndof)=0.D0
      endif   
 
      if (option%iblkfmt == 1) then
        blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
      endif
      p1=(natural_id_up)*option%ndof;p2=(natural_id_dn)*option%ndof
      do ii=0,option%ndof-1
        do jj=0,option%ndof-1
          if (local_id_up>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
            else
              blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
            endif
          endif
          if (local_id_dn>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p2+ii,p1+jj,-ra(ii+1,jj+1),ADD_VALUES,ierr)
            else
              blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
            endif
          endif
        enddo

        do jj=option%ndof,2*option%ndof-1
          if (local_id_up>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p1+ii,p2+jj-option%ndof,ra(ii+1,jj+1), &
                               ADD_VALUES,ierr)
            else
              blkmat12(ii+1,jj-option%ndof+1) = blkmat12(ii+1,jj-option%ndof+1) + &
                                                               ra(ii+1,jj+1)
            endif
          endif
          if (local_id_dn>0) then
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
        if (local_id_up>0) call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_up,blkmat11,ADD_VALUES,ierr)
        if (local_id_dn>0) call MatSetValuesBlocked(A,1,natural_id_dn,1,natural_id_dn,blkmat22,ADD_VALUES,ierr)
        if (local_id_up>0) call MatSetValuesBlocked(A,1,natural_id_up,1,natural_id_dn,blkmat12,ADD_VALUES,ierr)
        if (local_id_dn>0) call MatSetValuesBlocked(A,1,natural_id_dn,1,natural_id_up,blkmat21,ADD_VALUES,ierr)
      endif

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
  use Field_module
  use Option_module
 
  implicit none
  
  type(solution_type) :: solution 

 
  integer :: ierr
  integer :: i, index_var_begin,index_var_end
  integer :: p1
  integer :: iicap, iiphase
  integer :: local_id, ghosted_id

  PetscScalar, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_loc_p(:),&
                          var_loc_p(:), icap_loc_p(:),iphase_loc_p(:),ithrm_loc_p(:)
  
  real*8 :: dif(1:solution%option%nphase),res(1:solution%option%ndof)
 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  grid => solution%grid
  option => solution%option
  field => solution%field
 
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%accum, accum_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  !print *,'Richardsinitaccum  Gotten pointers'
 
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)    
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    iicap=int(icap_loc_p(ghosted_id))
    iiphase = int(iphase_loc_p(ghosted_id))
    dif(1)= option%difaq

    call pri_var_trans_Richards_ninc(yy_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase,&
                                option%scale,option%nphase,option%nspec, &
                                iicap , dif, &
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                size_var_node+size_var_use),option%itable,ierr, &
                                option%pref)
  enddo

!---------------------------------------------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...

    ghosted_id = grid%nL2G(local_id)
    
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    p1 = 1 + (local_id-1)*option%ndof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin-1 + size_var_use
    i = ithrm_loc_p(ghosted_id)
    
    call RichardsRes_ARCont(local_id,var_loc_p(index_var_begin:index_var_end), &
                            porosity_loc_p(ghosted_id),volume_p(local_id), &
                            option%dencpr(i),option,Res, &
                             0,ierr)
 

    accum_p(p1:p1+option%ndof-1)=Res(:) 
    
  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%accum, accum_p, ierr)
  call VecRestoreArrayF90(field%var_loc, var_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine pflow_Richards_initaccum


subroutine pflow_update_Richards(solution)

  use translator_Richards_module
  use pckr_module
  use Condition_module_old

  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
  use Coupler_module

  implicit none

  type(solution_type) :: solution 
    
  integer :: dof_offset
!geh added for transient boundary conditons
  integer :: nc, iithrm
  real*8 :: sw, pc(2), kr(2)
  integer :: ierr,iicap,iiphase, iiphase_old
  PetscScalar, pointer :: xx_p(:),icap_loc_p(:),ithrm_loc_p(:), &
                          iphase_loc_p(:), var_loc_p(:), yy_p(:), iphase_loc_old_p(:)
  real*8 :: dif(1:solution%option%nphase)
  integer :: local_id, ghosted_id        


  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  integer :: sum_connection  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => solution%grid
  option => solution%option
  field => solution%field
      
!geh added for transient boundary conditions      
  if (associated(field%imat) .and. option%iread_geom < 0) then
!commend out for now    call UpdateBoundaryConditions(option)
    yybc =field%xxbc
    vel_bc = field%velocitybc
  endif
!geh end
 
   ! if(ichange ==1)then
  call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%iphas_old_loc, iphase_loc_old_p, ierr)
  call VecGetArrayF90(field%var_loc,var_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    iicap = icap_loc_p(ghosted_id)
    iiphase = int(iphase_loc_p(ghosted_id))
    iiphase_old = int(iphase_loc_old_p(ghosted_id))
    dof_offset = (local_id-1)*option%ndof
    
     if(option%ndof>=3) then
        if (xx_p(dof_offset+3)<0.D0) xx_p(dof_offset+3)=1.D-6
     endif  
    
    dif(1) = option%difaq

    call pri_var_trans_Richards_ninc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap, dif,&
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                  size_var_node+size_var_use),&
                                option%itable,ierr, option%pref)
 
   enddo

  !geh added for transient boundary conditions  
  if (associated(field%imat) .and. option%iread_geom < 0) then

    boundary_condition => solution%boundary_conditions%first
    sum_connection = 0
    do 
      if (.not.associated(boundary_condition)) exit
    
      cur_connection_object => boundary_condition%connection

      do iconn = 1, cur_connection_object%num_connections
        sum_connection = sum_connection + 1
        nc = sum_connection

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
            sw= field%xxbc(1,nc)
            call pflow_pckr_richards_fw(iicap ,sw,pc,kr)    
            field%xxbc(1,nc) =  option%pref - pc(1)
          endif
      
          call pri_var_trans_Richards_ninc(field%xxbc(:,nc),field%iphasebc(nc), &
                                           option%scale,option%nphase,option%nspec, &
                                           iicap,dif, &
                                           field%varbc(1:size_var_use),option%itable,ierr, &
                                           option%pref)
        
          if (translator_check_cond_Richards(field%iphasebc(nc), &
                                             field%varbc(1:size_var_use), &
                                             option%nphase,option%nspec) /=1) then
            print *," Wrong bounday node init...  STOP!!!", field%xxbc(:,nc)
      
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
 
  call VecRestoreArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr);
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_old_loc, iphase_loc_old_p, ierr)
  call VecRestoreArrayF90(field%var_loc,var_loc_p,ierr)
   
  call translator_Richards_massbal(solution)
 ! endif 

  call VecCopy(field%xx, field%yy, ierr)   
  call VecCopy(field%iphas_loc, field%iphas_old_loc, ierr)   
   
  call  pflow_Richards_initaccum(solution)
! geh - comment
!translator_Richards_get_output is currently not necessary
!if uncommented, Vecs %pressure, %sat, %xmol must be created in pflow_init
!  call translator_Richards_get_output(grid,option)

end subroutine pflow_update_Richards





subroutine pflow_Richards_initadj(solution)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  use translator_Richards_module  
  use pckr_module, only: pflow_pckr_richards_fw

  use Connection_module
  use Solution_module
  use Grid_module
  use Option_module
  use Coupler_module
  
  implicit none

  type(solution_type) :: solution 

 
  integer :: ierr
  integer :: num_connection
  integer :: nc
  integer :: jn
  integer :: ii1,ii2,iicap
  integer :: iiphase,iithrm
  integer :: local_id, ghosted_id

  PetscScalar, pointer :: xx_p(:),var_loc_p(:)
  PetscScalar, pointer ::iphase_loc_p(:), ithrm_loc_p(:),icap_loc_p(:)
  
  real*8 :: dif(solution%option%nphase)
  real*8 :: pc(1:solution%option%nphase), kr(1:solution%option%nphase), sw

  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_object
  integer :: iconn
  integer :: sum_connection
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => solution%grid
  option => solution%option 
  field => solution%field 

  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%xx, xx_p, ierr)
  call VecGetArrayF90(field%var_loc, var_loc_p, ierr) 

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif

    jn = 1 + (local_id-1)*option%nphase
    ii1=1+(local_id-1)*option%nphase
    ii2=local_id*option%nphase
    iicap=int(icap_loc_p(local_id))
        
    iiphase = iphase_loc_p(local_id)
    dif(1)= option%difaq

   if(iiphase ==3)then

     sw= xx_p((local_id-1)*option%ndof+1)
     call pflow_pckr_richards_fw(iicap,sw,pc,kr)    
                 
     if(pc(1)>option%pcwmax(iicap))then
        print *,'INIT Warning: Pc>pcmax', sw, pc(1), iicap, option%pcwmax(iicap)
        pc(1)=option%pcwmax(iicap)
     endif 
     xx_p((local_id-1)*option%ndof+1)= option%pref - pc(1)
     print *,'Richards: Conv: ',local_id, sw, iicap, sw, pc(1),xx_p((local_id-1)*option%ndof+1:local_id*option%ndof)
    endif
    
    call pri_var_trans_Richards_ninc(xx_p((local_id-1)*option%ndof+1:local_id*option%ndof),iiphase, &
                                option%scale,option%nphase,option%nspec, &
                                iicap,  dif, &
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                  size_var_node+size_var_use), &
                                option%itable,ierr, option%pref)
   
    if (translator_check_cond_Richards(iiphase, &
                                        var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)* &
                                        size_var_node+size_var_use), &
                                        option%nphase,option%nspec) /= 1 ) then
      print *," Wrong internal node init...  STOP!!!"
      stop    
    endif 
  enddo

  boundary_condition => solution%boundary_conditions%first
  num_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit    
    num_connection = num_connection + boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo

  allocate(yybc(option%ndof,num_connection))
  allocate(vel_bc(option%nphase,num_connection))
  yybc =field%xxbc
  vel_bc = field%velocitybc

  boundary_condition => solution%boundary_conditions%first
  sum_connection = 0  
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_object => boundary_condition%connection
    
    do iconn = 1, cur_connection_object%num_connections
      sum_connection = sum_connection + 1
      nc = sum_connection
      
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
        iicap=int(icap_loc_p(local_id))
        iithrm=int(ithrm_loc_p(local_id)) 
        dif(1)= option%difaq
        
        if(field%iphasebc(nc) ==3)then
          sw= field%xxbc(1,nc)
          call pflow_pckr_richards_fw(iicap,sw,pc,kr)    
          field%xxbc(1,nc) =  option%pref - pc(1)
        endif

        
        
        call pri_var_trans_Richards_ninc(field%xxbc(:,nc),field%iphasebc(nc), &
                                         option%scale,option%nphase,option%nspec, &
                                         iicap, dif, &
                                         field%varbc(1:size_var_use),option%itable,ierr, &
                                         option%pref)
        
        if (translator_check_cond_Richards(field%iphasebc(nc), &
                                            field%varbc(1:size_var_use), &
                                            option%nphase,option%nspec) /=1) then
          print *," Wrong bounday node init...  STOP!!!", field%xxbc(:,nc)
        
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

  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%xx, xx_p, ierr)
  call VecRestoreArrayF90(field%var_loc, var_loc_p, ierr)
   
end subroutine pflow_Richards_initadj


subroutine createRichardsZeroArray(solution)

  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none

  type(solution_type) :: solution
  integer :: ncount, idof
  integer :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
    
  grid => solution%grid
  option => solution%option
  field => solution%field
  
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

end subroutine createRichardsZeroArray


end module Richards_module
