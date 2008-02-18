! introduced grid variables: e_total :: 1 dof
! translator_module                           c_total :: option%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (option%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

               
module MPHASE_module

  use mphase_field_module  ! located at top of translator_mixed_fluid_mph.F90
  use mphase_option_module ! located at top of translator_mixed_fluid_mph.F90

  implicit none
  
  private

#include "definitions.h"
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
!  PetscReal, parameter :: formeps   = 5.D-5
  PetscReal, parameter :: eps       = 1.D-5
  PetscReal, parameter :: floweps   = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0!1D-8
  PetscReal, parameter :: dfac = 1.D-8
!  PetscReal, parameter :: dfac = 1.D-2

  PetscInt,save :: size_var_use 
  PetscInt,save :: size_var_node
  PetscReal, allocatable,save :: Resold_AR(:,:), Resold_FL(:,:)
! Contributions to residual from accumlation/source/Reaction, flux(include diffusion)

  public MPHASEResidual, MPHASEJacobian, pflow_mphase_initaccum, &
         pflow_update_mphase,pflow_mphase_initadj, pflow_mphase_timecut,&
         pflow_mphase_setupini, MPhase_Update_Reason, MphaseSetup, &
         MphaseCheckpointWrite, MphaseCheckpointRead, MphaseMaxChange, &
         MphaseInitializeSolidReaction, MphaseGetVarFromArray, &
         MphaseGetTecplotHeader

  PetscInt, save :: n_zero_rows = 0
  PetscInt, pointer, save :: zero_rows_local(:)  ! 1-based indexing
  PetscInt, pointer, save :: zero_rows_local_ghosted(:) ! 0-based indexing
  
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
 
  PetscReal, pointer :: xx_p(:),yy_p(:)!,var_loc_p(:),iphase_loc_p(:)
  PetscInt :: re
  PetscErrorCode :: ierr
  PetscInt :: dof_offset
  PetscInt :: local_id
  !PetscInt :: re0, ierr, index, iipha
  !PetscReal, pointer :: sat(:),xmol(:)

  grid => realization%grid
  option => realization%option
  field => realization%field  

  mphase_option%iphch=0
    
  call VecGetArrayF90(field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr)
 !call VecGetArrayF90(field%var_loc, var_loc_p, ierr); 
 !call VecGetArrayF90(field%iphas, iphase_loc_p, ierr); 

  do local_id=1, grid%nlmax
    dof_offset=(local_id-1)*option%nflowdof
    do re = 1, option%nflowdof
!     if (option%snes_reason == -5) then 
!       xx_p(n0+re)= yy_p(n0+re)       !.5D0 * xx_p(n0+re) +.5D0 *yy_p(n0+re)
!     else
        xx_p(dof_offset+re)= yy_p(dof_offset+re)
!     endif
    enddo
  enddo 
  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr) 
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  
  !call VecCopy(field%flow_xx,field%flow_yy,ierr)
  !call pflow_mphase_initaccum(grid)
 
end subroutine pflow_mphase_timecut

! ************************************************************************** !
!
! MphaseCheckpointWrite: Writes vecs to checkpoint file
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseCheckpointWrite(grid, viewer)

  use Grid_module

  implicit none
  
  type(grid_type) :: grid
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call GridCreateVector(grid,VARDOF,global_var,GLOBAL)
  call GridLocalToGlobal(grid,mphase_field%var_loc,global_var,VARDOF)
  call VecView(global_var,viewer,ierr)
  call VecDestroy(global_var,ierr)
  
  ! solid volume fraction
  if (mphase_option%rk > 0.d0) then
    call VecView(mphase_field%phis, viewer, ierr)
  endif  
  
end subroutine MphaseCheckpointWrite

! ************************************************************************** !
!
! MphaseCheckpointRead: Reads vecs from checkpoint file
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseCheckpointRead(grid,viewer)

  use Grid_module

  implicit none
  
  type(grid_type) :: grid
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call GridCreateVector(grid,VARDOF,global_var,GLOBAL)
  call VecLoadIntoVector(viewer, global_var, ierr)
  call GridGlobalToLocal(grid,global_var,mphase_field%var_loc,VARDOF)
  call VecDestroy(global_var,ierr)
  ! solid volume fraction
  if (mphase_option%rk > 0.d0) then
    call VecLoadIntoVector(viewer, mphase_field%phis, ierr)
  endif  
  
end subroutine MphaseCheckpointRead
  
! ************************************************************************** !
!
! MphaseMaxChange: 
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseMaxChange(realization)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(grid_type), pointer :: grid
  
  PetscReal, pointer :: xx_p(:), yy_p(:), iphase_loc_p(:),var_loc_p(:),iphase_old_loc_p(:)
  PetscReal :: comp1,comp,cmp  
! PetscReal :: dsm,dcm  
  PetscReal :: dsm0,dcm0  
  PetscInt :: local_id, ghosted_id, dof_offset
! PetscInt :: j
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  grid => realization%grid
  
   call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
    call VecStrideNorm(field%flow_dxx,0,NORM_INFINITY,option%dpmax,ierr)
    call VecStrideNorm(field%flow_dxx,1,NORM_INFINITY,option%dtmpmax,ierr)

  call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  call VecGetArrayF90(field%iphas_old_loc, iphase_old_loc_p,ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  
  comp=0.D0;comp1=0.D0
  do local_id=1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    dof_offset=(local_id-1)*option%nflowdof 
    if(int(iphase_loc_p(ghosted_id)) == int(iphase_old_loc_p(ghosted_id)))then
     cmp=dabs(xx_p(dof_offset+3)-yy_p(dof_offset+3))
     if(int(iphase_loc_p(ghosted_id))==1 .or.int(iphase_loc_p(ghosted_id))==2)then
       if(comp<cmp) comp=cmp
     
    endif   
       if(int(iphase_loc_p(ghosted_id))==3)then
       if(comp1<cmp) comp1=cmp
      
    endif   
    else
!  print *,'phase changed', n, iphase_loc_p(n), iphase_old_p(n)

   endif
  enddo
  !call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_old_loc, iphase_old_loc_p,ierr)
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr)
 
  
  if(option%commsize >1)then
    call MPI_ALLREDUCE(comp1, dsm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dsm0,ONE_INTEGER, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(comp, dcm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dcm0,ONE_INTEGER, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    comp1 = dsm0
    comp = dcm0
  endif 

   option%dsmax=comp1
   option%dcmax=comp
!   print *, 'max change',grid%dpmax,option%dtmpmax,grid%dsmax,grid%dcmax

end subroutine MphaseMaxChange
        
! ************************************************************************** !
!
! MphaseSetup: Creates arrays for auxilliary variables, initializes, etc.
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseSetup(realization)

  use Realization_module
 
  implicit none
  
  type(realization_type) :: realization

  call MphaseInitializeSolidReaction(realization)

#if 0

! from option%iphch
    PetscInt :: iphch

  ! declarations formerly in field.F90
  PetscReal, pointer :: xphi_co2(:),xxphi_co2(:),den_co2(:), dden_co2(:)
  PetscReal, pointer :: xphi_co2_bc(:), xxphi_co2_bc(:)
  
  
  ! not sure if these are need in mphase...?
    Vec :: conc
    Vec :: ttemp, ttemp_loc, temp ! 1 dof  
    
    
   Vec :: var_loc

! someone will have to sort out whether these belong in mphase, from field.F90
    Vec :: ppressure, ppressure_loc, pressure, dp
    Vec :: ssat, ssat_loc, sat    ! saturation
    Vec :: xxmol, xxmol_loc, xmol ! mole fraction
    Vec :: density       ! Density at time k
    Vec :: ddensity, ddensity_loc  ! Density at time k+1
    Vec :: d_p, d_p_loc  ! dD/dp at time k+1
    Vec :: d_t, d_t_loc  ! dD/dT at time k+1
    Vec :: d_c, d_c_loc  ! dD/dT at time k+1
    Vec :: d_s, d_s_loc  ! dD/dT at time k+1
    Vec :: avgmw,avgmw_loc  ! Density at time k+1molecular weight at time k+1
    Vec :: avgmw_c,avgmw_c_loc
    Vec :: h             ! H     at time k
    Vec :: hh, hh_loc    ! H     at time k+1
    Vec :: h_p, h_p_loc  ! dH/dp at time k+1
    Vec :: h_t, h_t_loc  ! dH/dT at time k+1
    Vec :: h_c, h_c_loc  ! dD/dT at time k+1
    Vec :: h_s, h_s_loc  ! dD/dT at time k+1
    Vec :: u            ! H     at time k
    Vec :: uu, uu_loc    ! H     at time k+1
    Vec :: u_p, u_p_loc  ! dH/dp at time k+1
    Vec :: u_t, u_t_loc  ! dH/dT at time k+1
    Vec :: u_c, u_c_loc  ! dD/dT at time k+1
    Vec :: u_s, u_s_loc  ! dD/dT at time k+1
    Vec :: hen, hen_loc    ! H     at time k+1
    Vec :: hen_p, hen_p_loc  ! dH/dp at time k+1
    Vec :: hen_t, hen_t_loc  ! dH/dT at time k+1
    Vec :: hen_c, hen_c_loc  ! dD/dT at time k+1
    Vec :: hen_s, hen_s_loc  ! dD/dT at time k+1
    Vec :: df, df_loc    ! H     at time k+1
    Vec :: df_p, df_p_loc  ! dH/dp at time k+1
    Vec :: df_t, df_t_loc  ! dH/dT at time k+1
    Vec :: df_c, df_c_loc  ! dD/dT at time k+1
    Vec :: df_s, df_s_loc  ! dD/dT at time k+1
    Vec :: viscosity, viscosity_loc  !kept for early routine
   
    Vec :: v_p, v_p_loc  ! dv/dp at time k+1
    Vec :: v_t, v_t_loc  ! dv/dT at time k+1
    Vec :: pcw, pcw_loc    ! H     at time k+1
    Vec :: pc_p, pc_p_loc  ! dH/dp at time k+1
    Vec :: pc_t, pc_t_loc  ! dH/dT at time k+1
    Vec :: pc_c, pc_c_loc  ! dD/dT at time k+1
    Vec :: pc_s, pc_s_loc  ! dD/dT at time k+1
    Vec :: kvr, kvr_loc    ! H     at time k+1
    Vec :: kvr_p, kvr_p_loc  ! d/dp at time k+1
    Vec :: kvr_t, kvr_t_loc  ! dm/dT at time k+1
    Vec :: kvr_c, kvr_c_loc  ! d/d at time k+1
    Vec :: kvr_s, kvr_s_loc  ! dD/dT at time k+1

    PetscReal, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
    PetscReal, pointer :: vvlbc(:), vvgbc(:)
    PetscReal, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)

! from pflow_init
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)  
      call VecDuplicate(field%porosity_loc, field%ttemp_loc, ierr)
  end select

  ! should these be moved to their respective modules
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      ! nphase degrees of freedom
      call GridCreateVector(grid,NPHASEDOF,field%pressure,GLOBAL)
      call VecDuplicate(field%pressure, field%sat, ierr)
      call VecDuplicate(field%pressure, field%xmol, ierr)
      call VecDuplicate(field%pressure, field%ppressure, ierr)
      call VecDuplicate(field%pressure, field%ssat, ierr)
      call VecDuplicate(field%pressure, field%dp, ierr)
      call VecDuplicate(field%pressure, field%density, ierr)
      call VecDuplicate(field%pressure, field%ddensity, ierr)
      call VecDuplicate(field%pressure, field%avgmw, ierr)
      call VecDuplicate(field%pressure, field%d_p, ierr)
      call VecDuplicate(field%pressure, field%d_t, ierr)
      call VecDuplicate(field%pressure, field%h, ierr)
      call VecDuplicate(field%pressure, field%hh, ierr)
      call VecDuplicate(field%pressure, field%h_p, ierr)
      call VecDuplicate(field%pressure, field%h_t, ierr)
      call VecDuplicate(field%pressure, field%viscosity, ierr)
      call VecDuplicate(field%pressure, field%v_p, ierr)
      call VecDuplicate(field%pressure, field%v_t, ierr)
     ! xmol may not be nphase DOF, need change later 
      call VecDuplicate(field%pressure, field%flow_xxmol, ierr)
  end select

  ! should these be moved to their respective modules?
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      call GridCreateVector(grid,NPHASEDOF, field%ppressure_loc, LOCAL)
      call VecDuplicate(field%ppressure_loc, field%ssat_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%flow_xxmol_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%ddensity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%avgmw_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%hh_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%viscosity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_t_loc, ierr)
  end select

  ! move to mph
  select case(option%iflowmode)
    case(MPH_MODE)
      call GridCreateVector(grid,VARDOF, field%var_loc,LOCAL)
  end select

  ! move to mph
  select case(option%iflowmode)
    case(MPH_MODE)
      call initialize_span_wagner(option%itable,option%myrank)  
  end select

  ! move to mph
  select case(option%iflowmode)
    case(MPH_MODE)
      allocate(field%xphi_co2(grid%nlmax))
      allocate(field%xxphi_co2(grid%nlmax))
      allocate(field%den_co2(grid%nlmax))
      allocate(field%dden_co2(grid%nlmax))
      field%xphi_co2 = 1.d0
      field%xxphi_co2 = 1.d0
      field%den_co2 = 1.d0
      field%dden_co2 = 1.d0
  end select
  
  ! move to mph
! Note: VecAssemblyBegin/End needed to run on the Mac - pcl (11/21/03)!
  if (field%conc /= 0) then
    call VecAssemblyBegin(field%conc,ierr)
    call VecAssemblyEnd(field%conc,ierr)
  endif
  if (field%xmol /= 0) then
    call VecAssemblyBegin(field%xmol,ierr)
    call VecAssemblyEnd(field%xmol,ierr)
  endif


    
! These data types need to be localized within this mode, not in the field type
  if (option%iflowmode /= RICHARDS_MODE .and. &
      option%iflowmode /= RICHARDS_LITE_MODE) then
    allocate(field%xphi_co2_bc(temp_int))
    allocate(field%xxphi_co2_bc(temp_int))
  endif

  ! no longer supported in option, needs to be localized to mph  
  !set specific phase indices
  option%jh2o = 1; option%jgas =1
  select case(option%nphase)
    case(2)
      option%jco2 = 2
      option%jgas =3 
    case(3)
      option%jco2 = 2
      option%jgas =3 
  end select   

! from option.F90 (i.e. option%nvar) : no longer supported in option, needs to be localized to this module  
    PetscInt :: nvar
    
    
    
  ! these vecs need to be stored within this module, not in field
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      call VecDuplicate(field%porosity0, field%conc, ierr)
      call VecDuplicate(field%porosity0, field%temp, ierr)
      call VecDuplicate(field%porosity0, field%ttemp, ierr)
  end select    
  
   from field
   Vec :: phis
    field%phis = 0
  
#endif

  mphase_option%iphch=0
  
  call createMphaseZeroArray(realization)
  call pflow_mphase_initadj(realization)
  call pflow_update_mphase(realization)
  
end subroutine MphaseSetup

subroutine pflow_mphase_setupini(realization)
 
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Region_module
  use Structured_Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: initial_condition
  
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscInt :: local_id, ghosted_id, ibegin, iend, icell, idof
  PetscErrorCode :: ierr
  
  grid => realization%grid
  option => realization%option
  field => realization%field
    
  size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
  size_var_node = (option%nflowdof + 1) * size_var_use
  
  allocate(Resold_AR(grid%nlmax,option%nflowdof))
  allocate(Resold_FL(ConnectionGetNumberInList(grid%internal_connection_list), &
                     option%nflowdof))
  allocate(mphase_option%delx(option%nflowdof,grid%ngmax))
  mphase_option%delx=0.D0
   
  call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  
  initial_condition => realization%initial_conditions%first
  
  do
  
    if (.not.associated(initial_condition)) exit
    
    do icell=1,initial_condition%region%num_cells
      local_id = initial_condition%region%cell_ids(icell)
      ghosted_id = grid%nL2G(local_id)
      ibegin = (local_id-1)*option%nflowdof
      do idof=1,option%nflowdof
        xx_p(ibegin+idof) = &
          initial_condition%condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
      enddo
      iphase_loc_p(ghosted_id)=initial_condition%condition%iphase
    enddo
  
    initial_condition => initial_condition%next
  
  enddo
              
  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr)
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
   
  PetscInt, intent(out):: reason
  PetscReal, pointer :: xx_p(:),var_loc_p(:),iphase_loc_p(:), yy_p(:) !,r_p(:)
  PetscInt :: dof_offset, temp_reason
  PetscInt :: iipha
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id  

! PetscReal, pointer :: sat(:),xmol(:)
! PetscReal rmax(option%nflowdof)

  grid => realization%grid
  option => realization%option
  field => realization%field
  
  reason=1
 ! call SNESComputeFunction(option%snes,field%flow_xx,option%r,ierr)
 ! do n=1,option%nflowdof
 !  call VecStrideNorm(option%r,n-1,NORM_INFINITY,rmax(n),ierr)
 ! enddo
  
 ! if(rmax(1)>1.D0 .or. rmax(2)>1.D0 .or. rmax(3)>5.D0)then
 !   re=0;print *, 'Rmax error: ',rmax
 ! endif
  
  if(reason>0)then
  call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr); 
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); 
  
  do local_id = 1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif  
     dof_offset=(local_id-1)* option%nflowdof
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
!     n0=(n-1)* option%nflowdof
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
    call VecRestoreArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
    call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr) 
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr) 
  endif
 ! print *,' update reason', option%myrank, re,n,grid%nlmax
  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  if(option%commsize >1)then
    temp_reason = reason
    call MPI_ALLREDUCE(temp_reason,reason,ONE_INTEGER, MPI_INTEGER,MPI_SUM, &
    PETSC_COMM_WORLD,ierr)
  !print *,' update reason re'
    !call MPI_BCAST(re0,ONE_INTEGER, MPI_INTEGER, ZERO_INTEGER,PETSC_COMM_WORLD,ierr)
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
  PetscInt :: node_no
  PetscInt, optional:: ireac
  PetscErrorCode :: ierr
  PetscReal, target :: var_node(1:size_var_use)
  PetscReal :: Res_AR(1:option%nflowdof) 
  PetscReal :: vol,por,rock_dencpr
     
  PetscReal, pointer :: temp, pre_ref   ! 1 dof
  PetscReal, pointer :: sat(:), density(:), amw(:), h(:), u(:), pc(:), kvr(:) ! nphase dof
  PetscReal, pointer :: xmol(:), diff(:)            ! nphase*nspec
  
  PetscInt :: ibase, m,np, iireac=1
  PetscReal :: pvol,mol(option%nspec),eng
  
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
!  if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!   H2O
!    mol(1)= mol(1) - option%dt * option%rtot(node_no,1)
!   CO2
!    mol(2)= mol(2) - option%dt * option%rtot(node_no,2)
!   should include related energy change here
!  endif
  Res_AR(1:option%nflowdof-1)=mol(:)
  Res_AR(option%nflowdof)=eng

  nullify(temp, pre_ref, sat, density, amw, h,u, pc,kvr,xmol,diff)
  
end subroutine  MPHASERes_ARCont


subroutine MPHASERes_FLCont(nconn_no,area, &
                             var_node1,por1,tor1,sir1,dd1,perm1,Dk1,&
                             var_node2,por2,tor2,sir2,dd2,perm2,Dk2,dist_gravity,upweight,&
                             option, vv_darcy,Res_FL)
  
  use Option_module                             

  implicit none
  
  PetscInt :: nconn_no
  type(option_type) :: option
  PetscReal :: sir1(1:option%nphase),sir2(1:option%nphase)
  PetscReal, target :: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  PetscReal, target :: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  PetscReal :: por1,por2,tor1,tor2,perm1,perm2,Dk1,Dk2,dd1,dd2
  PetscReal :: vv_darcy(option%nphase),area
  PetscReal :: Res_FL(1:option%nflowdof) 
  
  PetscReal :: dist_gravity ! distance along gravity vector
     
  PetscReal, pointer :: temp1, pre_ref1   ! 1 dof
  PetscReal, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  PetscReal, pointer :: xmol1(:), diff1(:)            ! 
  
  PetscReal, pointer :: temp2, pre_ref2   ! 1 dof
  PetscReal, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  PetscReal, pointer :: xmol2(:), diff2(:)    
  
  PetscInt :: ibase, m,np, ind
  PetscReal :: fluxm(option%nspec),fluxe, v_darcy,q
  PetscReal :: uh,uxmol(1:option%nspec), ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond, gravity, dphi
  
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
      * dist_gravity
          
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
  Res_FL(1:option%nflowdof-1)=fluxm(:) * option%dt
  Res_FL(option%nflowdof)=fluxe * option%dt
 ! note: Res_FL is the flux contribution, for node 1 R = R + Res_FL
 ! 2 R = R - Res_FL
 !print *,'end FLcont'

  nullify(temp1, pre_ref1, sat1, density1, amw1, h1,u1, pc1,kvr1,xmol1,diff1)
  nullify(temp2, pre_ref2, sat2, density2, amw2, h2,u2, pc2,kvr2,xmol2,diff2)

end subroutine MPHASERes_FLCont


subroutine MPHASERes_FLBCCont(ibndtype,area,aux_vars, &
            var_node1,var_node2,por2,tor2,sir2,dd1,perm2,Dk2,dist_gravity,&
            option,field,vv_darcy,Res_FL)

  use Option_module
  use Field_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(option_type) :: option
  type(field_type) :: field
  PetscReal :: aux_vars(:)
  PetscReal :: dd1, sir2(1:option%nphase)
  PetscReal, target :: var_node1(1:2+7*option%nphase+2*option%nphase*option%nspec)
  PetscReal, target :: var_node2(1:2+7*option%nphase+2*option%nphase*option%nspec)
  PetscReal :: por2,perm2,Dk2,tor2
  PetscReal :: vv_darcy(option%nphase), area
  PetscReal :: Res_FL(1:option%nflowdof) 
  
  PetscReal :: dist_gravity   ! distance along gravity vector
     
  PetscReal, pointer :: temp1, pre_ref1   ! 1 dof
  PetscReal, pointer :: sat1(:), density1(:), amw1(:), h1(:), u1(:), pc1(:), kvr1(:)         ! nphase dof
  PetscReal, pointer :: xmol1(:), diff1(:)            ! 
  
  PetscReal, pointer :: temp2, pre_ref2   ! 1 dof
  PetscReal, pointer :: sat2(:), density2(:), amw2(:), h2(:), u2(:), pc2(:), kvr2(:)         ! nphase dof
  PetscReal, pointer :: xmol2(:), diff2(:)    
  
  PetscInt :: ibase,offset,iphase,idof,ispec,index
  PetscReal :: fluxm(option%nspec),fluxe,q(option%nphase)
  PetscReal :: uh(option%nphase),uxmol(option%nspec), ukvr,diff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave(option%nphase),cond,gravity, dphi

  PetscReal :: v_darcy
  
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

  select case (ibndtype(MPH_PRESSURE_DOF))
  
    case(DIRICHLET_BC)
    
      Dq= perm2 / dd1
      diffdp = por2*tor2/dd1*area
    ! Flow term
      do iphase =1,option%nphase
        if ((sat1(iphase) > sir2(iphase)) .or. (sat2(iphase) > sir2(iphase)))then
    
          upweight=1.D0
          if(sat1(iphase) <eps) then 
            upweight=0.d0
          else if(sat2(iphase) <eps) then 
            upweight=1.d0
          endif
          density_ave(iphase) = upweight*density1(iphase)+(1.D0-upweight)*density2(iphase)  
      
          gravity = (upweight*density1(iphase)*amw1(iphase) + &
          (1.D0-upweight)*density2(iphase)*amw2(iphase)) &
          * dist_gravity
          
          dphi = pre_ref1-pc1(iphase) - pre_ref2 + pc2(iphase) + gravity
    
!         if(option%iupstream(0,1) ==0)then
            if(dphi>=0.D0)then
              ukvr=kvr1(iphase)
              uh=h1(iphase)
              uxmol(1:option%nspec)=xmol1((iphase-1)*option%nspec+1 : iphase*option%nspec)
           else
              ukvr=kvr2(iphase)
              uh=h2(iphase)
              uxmol(1:option%nspec)=xmol2((iphase-1)*option%nspec+1 : iphase*option%nspec)
            endif      
!         else
!          if(option%iupstreambc(nbc_no, iphase) >=0 )then
!             ukvr=kvr1(iphase)
!             uh=h1(iphase)
!             uxmol(1:option%nspec)=xmol1((iphase-1)*option%nspec+1 : iphase*option%nspec)
!           else
!             ukvr=kvr2(iphase)
!             uh=h2(iphase)
!             uxmol(1:option%nspec)=xmol2((iphase-1)*option%nspec+1 : iphase*option%nspec)
!           endif      
!         endif  


          if(ukvr>floweps)then
            v_darcy= Dq * ukvr * dphi
            !option%vvl_loc(nbc_no) = v_darcy
            vv_darcy(iphase)=v_darcy
     
            q(iphase)=v_darcy * area
          
            do idof=1, option%nspec 
              fluxm(idof)=fluxm(idof) + q(iphase)*density_ave(iphase)*uxmol(idof)
            enddo
            fluxe = fluxe + q(iphase)*density_ave(iphase)*uh(iphase)
          endif
        endif
        
!       Diffusion term   
!       Note : average rule may not be correct  
        if ((sat1(iphase) > eps) .and. (sat2(iphase) > eps))then
          diff =diffdp * 0.25D0*(sat1(iphase)+sat2(iphase))*(density1(iphase)+density2(iphase))
          do idof = 1, option%nspec
            offset=idof+(iphase-1)*option%nspec
            fluxm(idof) = fluxm(idof) + diff * diff2(offset)&
            *( xmol1(offset)-xmol2(offset))
          enddo  
        endif
      enddo
      
!     conduction term

      Dk =  Dk2 / dd1
      cond=Dk*area*(temp1-temp2) 
      fluxe=fluxe + cond
 
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%nflowdof)=fluxe * option%dt

    case(NEUMANN_BC)
    
      if((dabs(aux_vars(MPH_PRESSURE_DOF))+ &
      dabs(aux_vars(option%nflowdof+MPH_PRESSURE_DOF)))>floweps)then
!       print *, 'FlowBC :',nbc_no,field%velocitybc(1,nbc_no),field%velocitybc(2,nbc_no)
        do iphase=1,option%nphase
 !        fluxm=0.D0; fluxe=0.D0
          v_darcy = aux_vars((iphase-1)*option%nflowdof+MPH_PRESSURE_DOF)
          vv_darcy(iphase) = aux_vars((iphase-1)*option%nflowdof+MPH_PRESSURE_DOF)
!         option%vvbc(j+(nc-2)*option%nphase)= field%velocitybc(j,nc)
      !   note different from 2 phase version

          if(v_darcy >0.d0)then 
            q(iphase) = v_darcy * density1(iphase) * area
            !q = 0.d0
            !flux = flux - q
            fluxe = fluxe - q(iphase)  * h1(iphase) 
            do idof=1, option%nspec
              fluxm(idof)=fluxm(idof) + q(iphase) * xmol1(idof + (iphase-1)*option%nspec)
            enddo 
          else 
            q(iphase) =  v_darcy * density2(iphase) * area   
            fluxe = fluxe - q(iphase)  * h2(iphase) 
            do idof=1, option%nspec
              fluxm(idof)=fluxm(idof) + q(iphase) * xmol2(idof + (iphase-1)*option%nspec)
            enddo 
          endif 
        enddo
      endif
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%nflowdof)=fluxe * option%dt

    case(HYDROSTATIC_BC)
   
      Dq= perm2 / dd1
      diffdp = por2*tor2/dd1*area
      
    ! Flow term
      do iphase =1,option%nphase
        if ((sat1(iphase) > sir2(iphase)) .or. (sat2(iphase) > sir2(iphase)))then
          upweight=1.D0
          if(sat1(iphase) <eps) then 
            upweight=0.d0
          else if(sat2(iphase) <eps) then 
            upweight=1.d0
          endif
          density_ave(iphase) = upweight*density1(iphase)+(1.D0-upweight)*density2(iphase)  
    
          gravity = (upweight*density1(iphase)*amw1(iphase) + &
          (1.D0-upweight)*density2(iphase)*amw2(iphase)) &
              * dist_gravity
        
          dphi = pre_ref1-pc1(iphase) - pre_ref2 + pc2(iphase) + gravity
    
!         if(option%iupstream(0,1) ==0)then
             if(dphi>=0.D0)then
               ukvr=kvr1(iphase)
               uh(iphase)=h1(iphase)
               uxmol(1:option%nspec)=xmol1((iphase-1)*option%nspec+1 : iphase*option%nspec)
            else
               ukvr=kvr2(iphase)
               uh(iphase)=h2(iphase)
               uxmol(1:option%nspec)=xmol2((iphase-1)*option%nspec+1 : iphase*option%nspec)
             endif      
!          else
!           if(option%iupstreambc(nbc_no, iphase) >=0 )then
!              ukvr=kvr1(iphase)
!              uh=h1(iphase)
!              uxmol(1:option%nspec)=xmol1((iphase-1)*option%nspec+1 : iphase*option%nspec)
!            else
!              ukvr=kvr2(iphase)
!              uh=h2(iphase)
!              uxmol(1:option%nspec)=xmol2((iphase-1)*option%nspec+1 : iphase*option%nspec)
!            endif      
!          endif  
 
         
          if(ukvr>floweps)then
            v_darcy= Dq * ukvr * dphi
               !option%vvl_loc(nbc_no) = v_darcy
            vv_darcy(iphase)=v_darcy
           
            q(iphase)=v_darcy * area
                
            do idof=1, option%nspec 
              fluxm(idof)=fluxm(idof) + q(iphase)*density_ave(iphase)*uxmol(idof)
            enddo
            fluxe = fluxe + q(iphase)*density_ave(iphase)*uh(iphase) 
          endif
        endif 
      enddo
   
      Res_FL(1:option%nspec)=fluxm(:)* option%dt
      Res_FL(option%nflowdof)=fluxe * option%dt
    
    case(4)
          
      Dk =  Dk2 / dd1
      cond = Dk*area*(temp1-temp2) 
      fluxe=fluxe + cond
     
      Res_FL(1:option%nspec)= 0.D0
      Res_FL(option%nflowdof)=fluxe * option%dt

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

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscInt :: ierr0
  PetscErrorCode :: ierr
  PetscInt :: nr
  PetscInt :: i, jn
  PetscInt :: ip2, p1, p2
  PetscInt :: ithrm1, ithrm2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  
  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               var_loc_p(:) 

  PetscReal :: xxbc(realization%option%nflowdof),varbc(1:size_var_use)
  PetscInt :: iphasebc, idof                    

  PetscReal, pointer :: iphase_loc_p(:),icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iicap,iiphase, index_var_begin, index_var_end,iicap1,iicap2,np

  PetscReal :: dd1, dd2, eng, &
!           eengl,eengg, &
!           fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl,&
!           fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q,  &
!           v_darcy,hflx,
            pvoldt, voldt, accum, pvol
  PetscReal :: dd, f1, f2, ff, perm1, perm2
! PetscReal :: Dphi,D0, por1, por2,density_ave
! PetscReal :: Dq, Dk  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
! PetscReal :: sat_pressure  ! Saturation pressure of water.
  PetscReal :: dw_kg, dw_mol,dif(realization%option%nphase)
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1!, qqsrc
! PetscReal :: cw1,cw2, xxlw,xxla,xxgw,xxga
  PetscReal :: cw
! PetscReal :: upweight
! PetscReal :: ukvr,uhh,uconc
  PetscReal :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp, rho, visc
  PetscReal :: Res(realization%option%nflowdof), vv_darcy(realization%option%nphase)
 
! PetscReal :: cond, den,
  PetscViewer :: viewer

  ! Connection variables
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity, upweight
  PetscInt :: ichange

  grid => realization%grid
  option => realization%option
  field => realization%field
 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  
  ierr = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
      
!   insure zero liquid sat not passed to ptran (no effect on pflow)
    if(xx_p((local_id-1)*option%nflowdof+3) < 0.D0)xx_p((local_id-1)*option%nflowdof+3) = zerocut
    if(xx_p((local_id-1)*option%nflowdof+3) > 1.D0)xx_p((local_id-1)*option%nflowdof+3) = 1.D0 - zerocut
    
!   check if p,T within range of table  
    if(xx_p((local_id-1)*option%nflowdof+1)< p0_tab*1D6 &
      .or. xx_p((local_id-1)*option%nflowdof+1)>(ntab_p*dp_tab + p0_tab)*1D6)then
       ierr=-1; exit  
     endif  
    if(xx_p((local_id-1)*option%nflowdof+2)< t0_tab -273.15D0 &
      .or. xx_p((local_id-1)*option%nflowdof+2)>ntab_t*dt_tab + t0_tab-273.15D0)then
      ierr=-1; exit
    endif
  enddo
  
  ierr0 = 0
  if(option%commsize >1)then
    call MPI_ALLREDUCE(ierr, ierr0,ONE_INTEGER, MPI_INTEGER,MPI_SUM, PETSC_COMM_WORLD,ierr)
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
  if(mphase_option%iphch>0 .and. mphase_option%iphch<=3)then
!  if(option%iphch<=3)then
    call Translator_MPhase_Switching(xx,realization,ZERO_INTEGER,ichange)   
  endif  
  mphase_option%iphch=mphase_option%iphch+1
   
  call VecRestoreArrayF90(xx, xx_p, ierr)

  call GridGlobalToLocal(grid,xx,field%flow_xx_loc,NFLOWDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)                          

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); CHKERRQ(ierr)

! there is potential possiblity that the pertubation of p may change the direction of pflow.
! once that happens, code may crash, namely wrong derive. 
  do ghosted_id = 1, grid%ngmax 

    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
      
    iiphase=int(iphase_loc_p(ghosted_id))

    mphase_option%delx(1,ghosted_id) = xx_loc_p((ghosted_id-1)*option%nflowdof+1)*dfac * 1.D-3
    mphase_option%delx(2,ghosted_id) = xx_loc_p((ghosted_id-1)*option%nflowdof+2)*dfac

    select case (iiphase)
      case (1)
       if(xx_loc_p((ghosted_id-1)*option%nflowdof+3) < .8)then
          mphase_option%delx(3,ghosted_id) =  dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3)
       else
          mphase_option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3) 
       endif
       if(mphase_option%delx(3,ghosted_id)<1D-9 .and. mphase_option%delx(3,ghosted_id)>=0.D0)mphase_option%delx(3,ghosted_id)=1D-9
       if(mphase_option%delx(3,ghosted_id)>-1D-9 .and. mphase_option%delx(3,ghosted_id)<0.D0)mphase_option%delx(3,ghosted_id)=-1D-9
      case(2)  
       if(xx_loc_p((ghosted_id-1)*option%nflowdof+3) <0.8)then
          mphase_option%delx(3,ghosted_id) =  dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3) 
       else
          mphase_option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3) 
       endif 
       if(mphase_option%delx(3,ghosted_id)<1D-9 .and. mphase_option%delx(3,ghosted_id)>=0.D0)mphase_option%delx(3,ghosted_id)=1D-9
       if(mphase_option%delx(3,ghosted_id)>-1D-9 .and. mphase_option%delx(3,ghosted_id)<0.D0)mphase_option%delx(3,ghosted_id)=-1D-9
      case(3)
        if(xx_loc_p((ghosted_id-1)*option%nflowdof+3) <=0.9)then
          mphase_option%delx(3,ghosted_id) = dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3) 
        else
          mphase_option%delx(3,ghosted_id) = -dfac*xx_loc_p((ghosted_id-1)*option%nflowdof+3) 
        endif 
        
       if(mphase_option%delx(3,ghosted_id)<1D-9 .and. mphase_option%delx(3,ghosted_id)>=0.D0)mphase_option%delx(3,ghosted_id)=1D-9
       if(mphase_option%delx(3,ghosted_id)>-1D-9 .and. mphase_option%delx(3,ghosted_id)<0.D0)mphase_option%delx(3,ghosted_id)=-1D-9
        
        if((mphase_option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%nflowdof+3))>1.D0)then
          mphase_option%delx(3,ghosted_id) = (1.D0-xx_loc_p((ghosted_id-1)*option%nflowdof+3))/1D5
        endif
        if((mphase_option%delx(3,ghosted_id)+xx_loc_p((ghosted_id-1)*option%nflowdof+3))<0.D0)then
          mphase_option%delx(3,ghosted_id) = xx_loc_p((ghosted_id-1)*option%nflowdof+3)/1D5
        endif
    end select
  enddo
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  
  
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
 ! print *, n,iicap,xx_p((n-1)*option%nflowdof+1:n*option%nflowdof) 
  !*******************************************
    call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof),iiphase,&
        option%scale,option%nphase,option%nspec,iicap, dif,&
    var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),&
    option%itable,option%m_nacl,ierr,mphase_field%xxphi_co2(local_id), mphase_field%dden_co2(local_id))



    if (option%ideriv .eq. 1) then
      call pri_var_trans_mph_winc(xx_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof),&
        mphase_option%delx(1:option%nflowdof,ghosted_id), iiphase,&
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
  call VecRestoreArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  ! call VecRestoreArrayF90(field%iphas_loce,iphase_loc_p,ierr)
  

  call GridLocalToLocal(grid,mphase_field%var_loc,mphase_field%var_loc,VARDOF)
  call GridLocalToLocal(grid,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)


! End distribute data 
! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
! call VecGetArrayF90(field%flow_yy, yy_p, ierr)
 

  ! notice:: here we assume porosity is constant
 
  call VecGetArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  call VecGetArrayF90(field%flow_yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'


  if (mphase_option%rk > 0.d0) then
    call VecGetArrayF90(mphase_field%phis,phis_p,ierr)
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
    
    p1 = 1 + (local_id-1)*option%nflowdof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    
    pvol = volume_p(local_id)*porosity_loc_p(ghosted_id)
    voldt = volume_p(local_id) / option%dt
    pvoldt = porosity_loc_p(ghosted_id) * voldt
    iiphase = iphase_loc_p(ghosted_id)
    i = ithrm_loc_p(ghosted_id)

    accum = 0.d0
    call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin: index_var_end),&
    porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(i), option, Res, ONE_INTEGER,ierr)
   
    r_p(p1:p1+option%nflowdof-1) = r_p(p1:p1+option%nflowdof-1) + Res(1:option%nflowdof)
    Resold_AR(local_id,1:option%nflowdof)= Res(1:option%nflowdof) 
  end do

!************************************************************************
! add source/sink terms
  source_sink => realization%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%condition%num_sub_conditions > MPH_CONCENTRATION_DOF) then
      enthalpy_flag = .true.
    else
      enthalpy_flag = .false.
    endif

    qsrc1 = source_sink%condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%condition%enthalpy%dataset%cur_value(1)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (enthalpy_flag) then
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%dt   
      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+2), &
                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
              
        r_p((local_id-1)*option%nflowdof + mphase_option%jh2o) = r_p((local_id-1)*option%nflowdof + mphase_option%jh2o) &
                                               - qsrc1 *option%dt
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - qsrc1*enth_src_h2o*option%dt
        Resold_AR(local_id,mphase_option%jh2o)= Resold_AR(local_id,mphase_option%jh2o) - qsrc1*option%dt
        Resold_AR(local_id,option%nflowdof)= Resold_AR(local_id,option%nflowdof) - qsrc1 * &
                                                             enth_src_h2o * option%dt
      endif  
    
      if (csrc1 > 0.d0) then ! injection

!        duan eos
!        call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!        call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!        enth_src_co2)
!        enth_src_co2=enth_src_co2 * 1.D-3     
 
        ! span-wagner
        rho = var_loc_p((ghosted_id-1)*size_var_node+4+option%nphase)*option%fmwco2 
        call co2_span_wagner(var_loc_p((ghosted_id-1)*size_var_node+2)*1.D-6,&
                             tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                             eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)

        ! units: rho [kg/m^3]; csrc1 [kmol/s]

        enth_src_co2 = enth_src_co2 * option%fmwco2
           
        r_p((local_id-1)*option%nflowdof + mphase_option%jco2) = r_p((local_id-1)*option%nflowdof + mphase_option%jco2) - csrc1*option%dt
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - csrc1 * enth_src_co2 *option%dt
        Resold_AR(local_id,mphase_option%jco2)= Resold_AR(local_id,mphase_option%jco2) - csrc1*option%dt
        Resold_AR(local_id,option%nflowdof)= Resold_AR(local_id,option%nflowdof) - csrc1 * enth_src_co2*option%dt

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
! option%iupstream = 0; option%iupstreambc = 0

  ! loop over internal connections
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or.  &
            field%imat(ghosted_id_dn) <= 0) cycle
      endif

      p1 = 1 + (local_id_up-1)*option%nflowdof 
      p2 = 1 + (local_id_dn-1)*option%nflowdof

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd1 = distance*fraction_upwind
      dd2 = distance-dd1 ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd2/(dd1+dd2)
        
      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm2 = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm1 = ithrm_loc_p(ghosted_id_up)
      ithrm2 = ithrm_loc_p(ghosted_id_dn)
      iicap1=int(icap_loc_p(ghosted_id_up))
      iicap2=int(icap_loc_p(ghosted_id_dn))
   
      D1 = option%ckwet(ithrm1)
      D2 = option%ckwet(ithrm2)

      call MPHASERes_FLCont(iconn,cur_connection_set%area(iconn), &
                            var_loc_p((ghosted_id_up-1)*size_var_node+1: &
                                      (ghosted_id_up-1)*size_var_node+size_var_use),&
                            porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                            option%sir(1:option%nphase,iicap1),dd1,perm1,D1,&
                            var_loc_p((ghosted_id_dn-1)*size_var_node+1: &
                                      (ghosted_id_dn-1)*size_var_node+size_var_use),&
                            porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                            option%sir(1:option%nphase,iicap2),dd2,perm2,D2,&
                            distance_gravity,upweight,option,vv_darcy,Res)
      field%internal_velocities(1,iconn) = vv_darcy(1) ! liquid
      field%internal_velocities(2,iconn) = vv_darcy(2) ! gas
    
      Resold_FL(iconn,1:option%nflowdof) = Res(1:option%nflowdof) 
    
      if(local_id_up>0)then
        r_p(p1:p1+option%nflowdof-1) = r_p(p1:p1+option%nflowdof-1) + Res(1:option%nflowdof)
      endif
   
      if(local_id_dn>0)then
        r_p(p2:p2+option%nflowdof-1) = r_p(p2:p2+option%nflowdof-1) - Res(1:option%nflowdof)
      endif

    enddo
  
    cur_connection_set => cur_connection_set%next

  end do

 
!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', option%myrank, option%sgbc    

 
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if(ghosted_id<=0)then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      end if

      p1 = 1 + (local_id-1) * option%nflowdof

      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))

! in order to match the legacy code
#define MATCH_LEGACY
#ifdef MATCH_LEGACY
      select case(boundary_condition%condition%itype(MPH_PRESSURE_DOF))
        case(DIRICHLET_BC)
          xxbc(1:option%nflowdof) = boundary_condition%aux_real_var(1:option%nflowdof,iconn)
          iphasebc = boundary_condition%aux_int_var(MPH_PRESSURE_DOF,iconn)
        case(NEUMANN_BC)
          xxbc(1:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
          iphasebc=int(iphase_loc_p(ghosted_id)) 
          if (boundary_condition%aux_real_var(MPH_PRESSURE_DOF,iconn) > 1.d-20) then
            xxbc(2:option%nflowdof) = boundary_condition%aux_real_var(2:option%nflowdof,iconn)
          endif
        case(HYDROSTATIC_BC)                              
          xxbc(1) = boundary_condition%aux_real_var(MPH_PRESSURE_DOF,iconn)
          xxbc(2:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+2:ghosted_id*option%nflowdof)
          iphasebc=int(iphase_loc_p(ghosted_id))
        case(ZERO_GRADIENT_BC)
          xxbc(1:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
          iphasebc=int(iphase_loc_p(ghosted_id))
        case(4)
          xxbc(MPH_PRESSURE_DOF) = xx_loc_p((ghosted_id-1)*option%nflowdof+MPH_PRESSURE_DOF)
          xxbc(3:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+3:ghosted_id*option%nflowdof)
          iphasebc=int(iphase_loc_p(ghosted_id))
      end select
#else
      do idof=1,option%nflowdof
        select case(boundary_condition%condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      select case(boundary_condition%condition%itype(MPH_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select
#endif      

      iicap=int(icap_loc_p(ghosted_id))  
       
    !*****************
      dif(1)= option%difaq
      dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    !*******************************************

  
      call pri_var_trans_mph_ninc(xxbc,iphasebc, &
                                  option%scale,option%nphase,option%nspec,iicap, &
                                  dif,varbc,option%itable, &
                                  option%m_nacl,ierr,mphase_field%xxphi_co2_bc(sum_connection),cw)
     
      call MPHASERes_FLBCCont(boundary_condition%condition%itype, &
                              cur_connection_set%area(iconn), &
                              boundary_condition%aux_real_var(:,iconn), &
                              varbc(1:size_var_use), &
                              var_loc_p((ghosted_id-1)*size_var_node+1: &
                                        (ghosted_id-1)*size_var_node+size_var_use), &
                              porosity_loc_p(ghosted_id),tor_loc_p(ghosted_id), &
                              option%sir(1:option%nphase,iicap), &
                              cur_connection_set%dist(0,iconn),perm1,D2, &
                              distance_gravity,option,field,vv_darcy,Res)
                              
      field%boundary_velocities(1,iconn) = vv_darcy(1)  ! liquid
      field%boundary_velocities(2,iconn) = vv_darcy(2)  ! gas

      r_p(p1:p1-1+option%nflowdof)= r_p(p1:p1-1+option%nflowdof) - Res(1:option%nflowdof)
      ResOld_AR(local_id,1:option%nflowdof) = ResOld_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
    enddo
    boundary_condition => boundary_condition%next
  enddo
 

!adjust residual to R/dt

  select case (mphase_option%idt_switch) 
    case(1) 
      r_p(:) = r_p(:)/option%dt
    case(-1)
      if(option%dt>1.D0) r_p(:) = r_p(:)/option%dt
  end select    
  
  do local_id = 1, grid%nlmax
    p1 = 1 + (local_id-1)*option%nflowdof
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
      p1 = 3 + (local_id-1)*option%nflowdof
      r_p(p1)=xx_loc_p(2 + (ghosted_id-1)*option%nflowdof)-yy_p(p1-1)
    enddo
  endif

  if (n_zero_rows > 0) then
    do i=1,n_zero_rows
      r_p(zero_rows_local(i)) = 0.
    enddo
  endif
   
!print *,'res closing pointer'
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  if (mphase_option%rk > 0.d0) then
    call VecRestoreArrayF90(mphase_field%phis,phis_p,ierr)
  endif

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'residual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'xx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

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
  use Field_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(inout) :: A, B
  type(realization_type) :: realization
 ! PetscInt, intent(inout) :: flag
  MatStructure flag

!   PetscInt :: j, jn, jm1, jm2,jmu,mu, 
  PetscInt :: i_upstream_revert
  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: i1, i2, jng, i
  PetscInt :: ithrm1, ithrm2
  PetscInt :: ip1, ip2 
  PetscInt :: p1,p2 
  PetscReal :: dum1, dum2
  PetscViewer :: viewer

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:),var_loc_p(:)
  PetscInt :: iicap,ii,jj,iiphas,iiphas1,iiphas2,iicap1,iicap2
  PetscInt :: index_var_begin, index_var_end
! PetscInt :: ibc_hencoeff
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho,dddt,dddp,fg,dfgdp,&
            dfgdt,eng,dhdt,dhdp,visc,dvdt,dvdp
! PetscReal :: cond, gravity, acc, density_ave, 
  PetscReal :: vv_darcy(realization%option%nphase), voldt, pvoldt
! PetscReal :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
!           flux, fluxh, fluxv, difff, diffg, diffl,
  PetscReal :: ff,dif(1:realization%option%nphase)
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd1, dd2, dd, f1, f2
! PetscReal :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
! PetscReal :: por1, por2, den
  PetscReal :: perm1, perm2
! PetscReal :: qu_rate, p_vapor,sat_pressure_t
! PetscReal :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
! PetscReal :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! PetscReal :: sat_pressure  ! Saturation pressure of water.
! PetscReal :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu, sat_ave
  PetscReal :: ra(1:realization%option%nflowdof,1:2*realization%option%nflowdof)  
! PetscReal :: uhh, uconc, ukvr
! PetscReal :: upweight,m1weight,m2weight,mbweight,mnweight
  PetscReal :: delxbc(1:realization%option%nflowdof)
  PetscReal :: blkmat11(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            blkmat12(1:realization%option%nflowdof,1:realization%option%nflowdof),&
            blkmat21(1:realization%option%nflowdof,1:realization%option%nflowdof),&
            blkmat22(1:realization%option%nflowdof,1:realization%option%nflowdof)
  PetscReal :: ResInc(1:realization%grid%nlmax,1:realization%option%nflowdof,1:realization%option%nflowdof), &
               res(1:realization%option%nflowdof)  
  PetscReal :: max_dev, norm
  PetscReal :: xxbc(realization%option%nflowdof), varbc(1:size_var_node)
  PetscInt :: iphasebc, idof
  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt ::  natural_id_up,natural_id_dn
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity, upweight
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

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
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
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr)

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
    do nvar=1, option%nflowdof
   
      index_var_begin=(ghosted_id-1)*size_var_node+nvar*size_var_use+1
      index_var_end = index_var_begin -1 + size_var_use

      call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin : index_var_end),&
        porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(int(ithrm_loc_p(ghosted_id))),&
        option, Res,ONE_INTEGER,ierr)
      
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
! add source/sink terms
  source_sink => realization%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit

! not needed when computing derivative since enthalpy is independent of solution vars    
    ! check whether enthalpy dof is included
!    if (size(source_sink%condition%cur_value) > MPH_CONCENTRATION_DOF) then
!      enthalpy_flag = .true.
!    else
!      enthalpy_flag = .false.
!    endif

    qsrc1 = source_sink%condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%condition%concentration%dataset%cur_value(1)
!    if (enthalpy_flag) hsrc1 = source_sink%condition%cur_value(MPH_ENTHALPY_DOF)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (qsrc1 > 0.d0) then ! injection
      
        do nvar=1,option%nflowdof      
          call wateos_noderiv(tsrc1,var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+2),&
                              dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)

!         units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!         qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
          ResInc(local_id,mphase_option%jh2o,nvar)=  ResInc(local_id,mphase_option%jh2o,nvar) - qsrc1*option%dt
          ResInc(local_id,option%nflowdof,nvar)=  ResInc(local_id,option%nflowdof,nvar) - qsrc1*enth_src_h2o*option%dt
        enddo

      endif  
    
      if (csrc1 > 0.d0) then ! injection

!       duan eos
!       call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!       call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!       enth_src_co2)
!       enth_src_co2=enth_src_co2 * 1.D-3     
 
        ! span-wagner
        do nvar=1,option%nflowdof     
          rho = var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+4+option%nphase)*option%fmwco2 
          call co2_span_wagner(var_loc_p((ghosted_id-1)*size_var_node+nvar*size_var_use+2)*1.D-6,&
                               tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                               eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)

          ! units: rho [kg/m^3]; csrc1 [kmol/s]

          enth_src_co2 = enth_src_co2 * option%fmwco2

          ResInc(local_id,mphase_option%jco2,nvar)=  ResInc(local_id,mphase_option%jco2,nvar) - csrc1*option%dt
          ResInc(local_id,option%nflowdof,nvar)=  ResInc(local_id,option%nflowdof,nvar) - csrc1*enth_src_co2*option%dt

        enddo

      endif
  
  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
      
  !  endif
    enddo
    source_sink => source_sink%next
  enddo


  ! print *,' Mph Jaco Finished source terms'
  
! Contribution from BC
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if(ghosted_id<=0)then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      end if
  
      p1 = 1 + (local_id-1) * option%nflowdof
      
      ithrm2 = ithrm_loc_p(ghosted_id)
      D2 = option%ckwet(ithrm2)

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
       
      delxbc=0.D0

! in order to match the legacy code
#ifdef MATCH_LEGACY
      select case(boundary_condition%condition%itype(MPH_PRESSURE_DOF))
        case(DIRICHLET_BC)
          xxbc(1:option%nflowdof) = boundary_condition%aux_real_var(1:option%nflowdof,iconn)
          delxbc(1:option%nflowdof) = 0.d0
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC)
          xxbc(1:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
          delxbc(1:option%nflowdof) = mphase_option%delx(1:option%nflowdof,ghosted_id) 
          iphasebc=int(iphase_loc_p(ghosted_id))                               
          if (boundary_condition%aux_real_var(MPH_PRESSURE_DOF,iconn) > 1.d-20) then
            xxbc(2:option%nflowdof) = boundary_condition%aux_real_var(2:option%nflowdof,iconn)
            delxbc(2:option%nflowdof) = 0.d0
          endif
        case(HYDROSTATIC_BC)
          xxbc(1) = boundary_condition%aux_real_var(MPH_PRESSURE_DOF,iconn)
          xxbc(2:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+2:ghosted_id*option%nflowdof)
          delxbc(1) = 0.d0
          delxbc(2:option%nflowdof) = mphase_option%delx(2:option%nflowdof,ghosted_id) 
          iphasebc=int(iphase_loc_p(ghosted_id))
        case(ZERO_GRADIENT_BC)
          xxbc(1:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
          delxbc(1:option%nflowdof) = mphase_option%delx(1:option%nflowdof,ghosted_id) 
          iphasebc=int(iphase_loc_p(ghosted_id))
        case(4)
          xxbc(1) = xx_loc_p((ghosted_id-1)*option%nflowdof+MPH_PRESSURE_DOF)
          xxbc(3:option%nflowdof) = xx_loc_p((ghosted_id-1)*option%nflowdof+3:ghosted_id*option%nflowdof)
          delxbc(1) = mphase_option%delx(1,ghosted_id) 
          delxbc(3:option%nflowdof) = mphase_option%delx(3:option%nflowdof,ghosted_id) 
          iphasebc=int(iphase_loc_p(ghosted_id))
      end select
#else
      do idof=1,option%nflowdof
        select case(boundary_condition%condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
            delxbc(idof) = 0.d0
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            delxbc(idof) = mphase_option%delx(idof,ghosted_id) 
        end select
      enddo
      
      select case(boundary_condition%condition%itype(MPH_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select
#endif      

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
      call pri_var_trans_mph_ninc(xxbc,iphasebc,&
                              option%scale,option%nphase,option%nspec,iicap,dif,&
                              varbc(1:size_var_use),option%itable, &
                              option%m_nacl,ierr, dum1, dum2)
  
      call pri_var_trans_mph_winc(xxbc,delxbc,&
                              iphasebc,option%scale,option%nphase, &
                              option%nspec,iicap,dif(1:option%nphase), &
                              varbc(size_var_use+1:(option%nflowdof+1)*size_var_use), &
                              option%itable,option%m_nacl,ierr)
            
!    print *,' Mph Jaco BC terms: finish increment'
      do nvar=1,option%nflowdof
   
        call MPHASERes_FLBCCont(boundary_condition%condition%itype, &
                                cur_connection_set%area(iconn), &
                                boundary_condition%aux_real_var(:,iconn), &
                                varbc(nvar*size_var_use+1: &
                                           (nvar+1)*size_var_use), &
                                var_loc_p((ghosted_id-1)*size_var_node+ &
                                             nvar*size_var_use+1: &
                                          (ghosted_id-1)*size_var_node+ &
                                             nvar*size_var_use+size_var_use), &
                                porosity_loc_p(ghosted_id),tor_loc_p(ghosted_id), &
                                option%sir(1:option%nphase,iicap), &
                                cur_connection_set%dist(0,iconn),&
                                perm1,D2,distance_gravity,option,field,vv_darcy,Res)

    
        ResInc(local_id,1:option%nflowdof,nvar) = ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
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
    p1 = (ghosted_id-1)*option%nflowdof ! = 1 + (ng-1)*option%nflowdof-1
   
    max_dev=0.D0
    do neq=1, option%nflowdof
      do nvar=1, option%nflowdof
        ra(neq,nvar)=ResInc(local_id,neq,nvar)/mphase_option%delx(nvar,ghosted_id) - &
                     ResOld_AR(local_id,neq)/mphase_option%delx(nvar,ghosted_id)
          if(max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
     
      enddo      
    enddo
    if(option%use_isoth==PETSC_TRUE)then
      ra(:,2)=0.D0
      ra(3,1:option%nflowdof)=0.D0
      ra(3,2)=1.D0
    endif     

   
  ! if(max_dev<1D-5)then
  !  print *,'Mph Jaco max dev = ', max_dev
 !  endif
  
    select case(mphase_option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%dt
      case(-1)
        if(option%dt>1) ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%dt
    end select
      

    if (option%iblkfmt == 0) then
      p1=(natural_id_up)*option%nflowdof
     
      do ii=0,option%nflowdof-1
        do jj=0,option%nflowdof-1
          call MatSetValue(A,p1+ii,p1+jj,ra(ii+1,jj+1)/ volume_p(local_id),ADD_VALUES,ierr)
        enddo
       enddo
    else
      !ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%dt
      blkmat11=ra(1:option%nflowdof,1:option%nflowdof)
    
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
  cur_connection_set => connection_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or. field%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      natural_id_up = grid%nG2N(ghosted_id_up)
      natural_id_dn = grid%nG2N(ghosted_id_dn)
      p2 =  (ghosted_id_dn-1)*option%nflowdof
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd1 = distance*fraction_upwind
      dd2 = distance-dd1 ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd2/(dd1+dd2)
    
      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm2 = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      iiphas1 = iphase_loc_p(ghosted_id_up)
      iiphas2 = iphase_loc_p(ghosted_id_dn)

      ithrm1 = ithrm_loc_p(ghosted_id_up)
      ithrm2 = ithrm_loc_p(ghosted_id_dn)
      D1 = option%ckwet(ithrm1)
      D2 = option%ckwet(ithrm2)
    
      iicap1 = int(icap_loc_p(ghosted_id_up))
      iicap2 = int(icap_loc_p(ghosted_id_dn))
  
  ! do neq = 1, option%nflowdof
      do nvar = 1, option%nflowdof
    
        call MPHASERes_FLCont(sum_connection,cur_connection_set%area(iconn), &
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

        ra(:,nvar)= Res(:)/mphase_option%delx(nvar,ghosted_id_up)-ResOld_FL(sum_connection,:)/mphase_option%delx(nvar,ghosted_id_up)

!     if(vv_darcy(1)>0.D0 .and. option%iupstream(sum_connection,1) == -1) i_upstream_revert =1 
!     if(vv_darcy(1)<0.D0 .and. option%iupstream(sum_connection,1) == 1)  i_upstream_revert =1
!     if(vv_darcy(2)>0.D0 .and. option%iupstream(sum_connection,2) == -1) i_upstream_revert =2
!     if(vv_darcy(2)<0.D0 .and. option%iupstream(sum_connection,2) == 1 ) i_upstream_revert =2
  
       
        call MPHASERes_FLCont(sum_connection,cur_connection_set%area(iconn), &
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
 
        ra(:,nvar+option%nflowdof)=Res(:)/mphase_option%delx(nvar,ghosted_id_dn) &
                              -ResOld_FL(sum_connection,:)/mphase_option%delx(nvar,ghosted_id_dn)

     
!     if(vv_darcy(1)>0.D0 .and. option%iupstream(sum_connection,1) == -1) i_upstream_revert =3 
!     if(vv_darcy(1)<0.D0 .and. option%iupstream(sum_connection,1) == 1)  i_upstream_revert =3
!     if(vv_darcy(2)>0.D0 .and. option%iupstream(sum_connection,2) == -1) i_upstream_revert =4
!     if(vv_darcy(2)<0.D0 .and. option%iupstream(sum_connection,2) == 1 ) i_upstream_revert =4
   
      enddo
  
   !   print *,' Mph Jaco Finished NC terms'
  
  ! enddo
   
      if (option%use_isoth==PETSC_TRUE) then
        ra(3,1:2*option%nflowdof)=0.D0
        ra(:,2)=0.D0
        ra(:,2+option%nflowdof)=0.D0
      endif   
 

      if (option%iblkfmt == 1) then
        blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
      endif
    
      p1=(natural_id_up)*option%nflowdof
      p2=(natural_id_dn)*option%nflowdof
 
      select case(mphase_option%idt_switch)
        case(1)
          ra =ra / option%dt
        case(-1)  
          if (option%dt>1) ra =ra / option%dt
      end select  
       
      do ii=0,option%nflowdof-1
        do jj=0,option%nflowdof-1
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
   
        do jj=option%nflowdof,2*option%nflowdof-1
          if (local_id_up>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p1+ii,p2+jj-option%nflowdof,ra(ii+1,jj+1)/volume_p(local_id_up),ADD_VALUES,ierr)
            else
              blkmat12(ii+1,jj-option%nflowdof+1) = blkmat12(ii+1,jj-option%nflowdof+1) + ra(ii+1,jj+1)
            endif
          endif
          if (local_id_dn>0) then
            if (option%iblkfmt == 0) then
              call MatSetValue(A,p2+ii,p2+jj-option%nflowdof,-ra(ii+1,jj+1)/volume_p(local_id_dn),ADD_VALUES,ierr)
            else
              blkmat22(ii+1,jj-option%nflowdof+1) =  blkmat22(ii+1,jj-option%nflowdof+1) - ra(ii+1,jj+1)
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
    cur_connection_set => cur_connection_set%next
  enddo
 

  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (mphase_option%rk > 0.d0) then
    call VecRestoreArrayF90(mphase_field%phis,phis_p,ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  f1 = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,f1,ierr) 
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    p1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      p2 = p1-1
    elseif (ii == option%nflowdof-1) then
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

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr)
    if (option%myrank == 0) print *, '1 norm:', norm
    call MatNorm(A,NORM_FROBENIUS,norm,ierr)
    if (option%myrank == 0) print *, '2 norm:', norm
    call MatNorm(A,NORM_INFINITY,norm,ierr)
    if (option%myrank == 0) print *, 'inf norm:', norm
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
!    if (option%myrank == 0) print *, 'max:', i, norm
  endif

end subroutine MPHASEJacobian



subroutine pflow_mphase_initaccum(realization)
 
  use translator_mph_module, only : pri_var_trans_mph_ninc
  
  use Realization_module
  use Grid_module
  use Field_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization 
  
  PetscErrorCode :: ierr
  PetscInt :: i, index_var_begin,index_var_end
  PetscInt :: p1
  PetscInt :: iicap, iiphase
  PetscInt :: local_id, ghosted_id  

  PetscReal, pointer :: accum_p(:),yy_p(:),volume_p(:),porosity_loc_p(:),&
                          var_loc_p(:), icap_loc_p(:),iphase_loc_p(:),ithrm_loc_p(:)
  
 !PetscInt, pointer ::iphase_loc_p(:)
  
! PetscReal :: sat_pressure
  PetscReal :: pvol, satw  ! Saturation pressure of water.
  PetscReal :: dif(1:realization%option%nphase),res(1:realization%option%nflowdof)
 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p,ierr)
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

    call pri_var_trans_mph_ninc(yy_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof),iiphase,&
                                option%scale,option%nphase,option%nspec, iicap, dif,&
                                var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),option%itable,&
                                option%m_nacl,ierr, satw, pvol)

  enddo

  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p,ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p,ierr)

!---------------------------------------------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif 
    p1 = 1 + (local_id-1)*option%nflowdof
    index_var_begin=(ghosted_id-1)*size_var_node+1
    index_var_end = index_var_begin -1 + size_var_use
    i = ithrm_loc_p(ghosted_id)
    
    call MPHASERes_ARCont(local_id, var_loc_p(index_var_begin: index_var_end),&
                          porosity_loc_p(ghosted_id),volume_p(local_id),option%dencpr(i), option, &
                          Res, ZERO_INTEGER,ierr)
 

    accum_p(p1:p1+option%nflowdof-1)=Res(:) 

   !print *, 'init m accum ', n,  Res 

! print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*option%nspec+(n-1)*option%nphase*option%nspec),kvr_p(jn),kvr_p(jn+1)

  enddo

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine pflow_mphase_initaccum


subroutine pflow_update_mphase(realization)
  
  use translator_mph_module, only : pri_var_trans_mph_ninc, &
                                    translator_mph_get_output, &
                                    translator_mphase_massbal, &
                                    translator_check_phase_cond
  
  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module 
  use Field_module 

  implicit none

  type(realization_type) :: realization 
    
  PetscInt :: dof_offset
  PetscInt :: iithrm
  PetscInt :: iicap,iiphase
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_p(:),icap_loc_p(:),ithrm_loc_p(:),iphase_loc_p(:), var_loc_p(:), phis_p(:)
  PetscReal :: dif(1:realization%option%nphase), dum1, dum2           
  PetscInt :: local_id, ghosted_id     
  PetscReal :: xxbc(realization%option%nflowdof), varbc(size_var_use)
  PetscInt :: iphasebc, idof, n
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => realization%grid
  option => realization%option
  field => realization%field
        
  ! if (option%rk > 0.d0) call Rock_Change(grid)
  ! call Translator_MPhase_Switching(field%flow_xx,grid,1,ierr)
  ! print *,'MPhase_Update done'
 
   !if(ichange ==1)then
    call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
    call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
    call VecGetArrayF90(mphase_field%var_loc,var_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)    
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
    
    iicap = icap_loc_p(ghosted_id)
    iiphase = iphase_loc_p(ghosted_id)
    dof_offset=(local_id-1)*option%nflowdof
    if(xx_p(dof_offset+3)<0.D0) xx_p(dof_offset+3)=zerocut

    ! *****************
      dif(1)= option%difaq
      dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
    ! *******************************************
      call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof),iiphase,&
                                  option%scale,option%nphase,option%nspec, iicap, dif,&
                                  var_loc_p((ghosted_id-1)*size_var_node+1:(ghosted_id-1)*size_var_node+size_var_use),&
                                  option%itable,option%m_nacl,ierr, dum1, dum2)

    enddo

  !geh added for transient boundary conditions  
  if (associated(field%imat) .and. option%iread_geom < 0) then

    boundary_condition => realization%boundary_conditions%first
    sum_connection = 0
    do 
      if (.not.associated(boundary_condition)) exit
    
      cur_connection_set => boundary_condition%connection

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        if (associated(field%imat)) then
          if (field%imat(ghosted_id) <= 0) cycle
        endif
       
        if (local_id<0) then
          print *, "Wrong boundary node index... STOP!!!"
          stop
        endif

        if (boundary_condition%condition%itype(MPH_PRESSURE_DOF) == DIRICHLET_BC) then
          do idof=1,option%nflowdof
            select case(boundary_condition%condition%itype(idof))
              case(DIRICHLET_BC,HYDROSTATIC_BC)
                xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
              case(NEUMANN_BC,ZERO_GRADIENT_BC)
                xxbc(idof) = xx_p((local_id-1)*option%nflowdof+idof)
            end select
          enddo
      
          select case(boundary_condition%condition%itype(MPH_PRESSURE_DOF))
            case(DIRICHLET_BC,HYDROSTATIC_BC)
              iphasebc = boundary_condition%aux_int_var(1,iconn)
            case(NEUMANN_BC,ZERO_GRADIENT_BC)
              iphasebc=int(iphase_loc_p(ghosted_id))                
          end select

          iicap = icap_loc_p(ghosted_id)
          dof_offset=(local_id-1)*option%nflowdof
          if(xx_p(dof_offset+3)<0.D0) xx_p(dof_offset+3)=zerocut

          ! *****************
          dif(1)= option%difaq
          dif(2)= option%cdiff(int(ithrm_loc_p(ghosted_id)))
          ! *****************
          call pri_var_trans_mph_ninc(xxbc,iphasebc,option%scale,option%nphase,option%nspec, &
                                      iicap,dif,&
                                      varbc,option%itable,option%m_nacl,ierr,dum1, dum2)

          if (translator_check_phase_cond(iphasebc,varbc,option%nphase,option%nspec) /=1) then
            print *," Wrong bounday node init...  STOP!!!", xxbc
      
            print *, varbc
            stop    
          endif 
        endif
      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif

   
    call VecRestoreArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)  
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
    call VecRestoreArrayF90(mphase_field%var_loc,var_loc_p,ierr)
     
    if(option%nphase>1) call translator_mphase_massbal(realization)
   ! endif 

    call VecCopy(field%flow_xx, field%flow_yy, ierr)   
    call VecCopy(field%iphas_loc, field%iphas_old_loc, ierr)   
     
    call  pflow_mphase_initaccum(realization)
      !print *,'pflow_mphase_initaccum done'
    call translator_mph_get_output(realization)
 !  print *,'translator_get_output done'
  ! the output variables should be put into option%pressure, temp,xmol,sat...
  ! otherwise need to rewrite the pflow_output

  !integrate solid volume fraction using explicit finite difference
  if (mphase_option%rk > 0.d0) then
    call VecGetArrayF90(mphase_field%phis,phis_p,ierr)
    do n = 1, grid%nlmax
      phis_p(n) = phis_p(n) + option%dt * mphase_option%vbars * mphase_option%rate(n)
      if (phis_p(n) < 0.d0) phis_p(n) = 0.d0
      mphase_option%area_var(n) = (phis_p(n)/mphase_option%phis0)**mphase_option%pwrsrf
      
!     print *,'update: ',n,phis_p(n),option%rate(n),grid%area_var(n)
    enddo
    call VecRestoreArrayF90(mphase_field%phis,phis_p,ierr)
  endif  

end subroutine pflow_update_mphase


subroutine pflow_mphase_initadj(realization)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  ! geh - will not compile without the 'only:' statement
  use translator_mph_module, only : pri_var_trans_mph_ninc, &
                                    translator_check_phase_cond

  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  
  implicit none
 
  type(realization_type) :: realization 

  PetscErrorCode :: ierr
  PetscInt :: num_connection  
  PetscInt :: jn
  PetscInt :: iicap
  PetscInt :: iiphase,iithrm
  PetscInt :: local_id, ghosted_id 

  PetscReal, pointer :: xx_p(:),var_loc_p(:)
  PetscReal, pointer ::iphase_loc_p(:), ithrm_loc_p(:),icap_loc_p(:)
  
  PetscReal  dif(realization%option%nphase), dum1, dum2
  
! PetscReal :: temp1
! PetscReal, parameter :: Rg=8.31415D0
  PetscReal :: xxbc(realization%option%nflowdof),varbc(1:size_var_use)
  PetscInt :: iphasebc, idof 

  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  grid => realization%grid
  option => realization%option 
  field => realization%field 


  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecGetArrayF90(field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr)
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
    call pri_var_trans_mph_ninc(xx_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof),iiphase,&
      option%scale,option%nphase,option%nspec, iicap,  dif,&
      var_loc_p((ghosted_id-1)*size_var_node+1: (ghosted_id-1)*size_var_node+size_var_use), &
      option%itable,option%m_nacl,ierr, dum1, dum2)
   
   !print *, xx_p((n-1)*option%nflowdof+1:n*option%nflowdof)
    if(translator_check_phase_cond(iiphase, &
      var_loc_p((ghosted_id-1)*size_var_node+1: (ghosted_id-1)*size_var_node+size_var_use),&
      option%nphase,option%nspec) /= 1 ) then
      print *," Wrong internal node init...  STOP!!!"
      stop    
    endif 
  enddo


! implemented in richards, but not here????
!  boundary_condition => realization%boundary_conditions%first
!  num_connection = 0
!  do 
!    if (.not.associated(boundary_condition)) exit    
!    num_connection = num_connection + boundary_condition%connection%num_connections
!    boundary_condition => boundary_condition%next
!  enddo

!  allocate(yybc(option%nflowdof,num_connection))
!  allocate(vel_bc(option%nphase,num_connection))
!  yybc =field%flow_xxbc
!  vel_bc = field%velocitybc

  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0  
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      
      local_id = cur_connection_set%id_dn(iconn)
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

        call pri_var_trans_mph_ninc(xxbc,iphasebc,&
                                    option%scale,option%nphase,option%nspec,iicap, &
                                    dif,varbc, &
                                    option%itable,option%m_nacl,ierr, dum1, dum2)
      
        if (translator_check_phase_cond(iphasebc, &
                                        varbc, &
                                        option%nphase,option%nspec) &
            /=1) then
          print *," Wrong bounday node init...  STOP!!!", xxbc
      
          print *, varbc
          stop    
        endif 
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  !print *,kgjkdf
  
  !call VecCopy(field%iphas,field%iphas_old,ierr)

end subroutine pflow_mphase_initadj



subroutine createMphaseZeroArray(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id

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
        n_zero_rows = n_zero_rows + option%nflowdof
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
        do idof = 1, option%nflowdof
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof+idof-1
        enddo
      else
#ifdef ISOTHERMAL
        ncount = ncount + 1
        zero_rows_local(ncount) = local_id*option%nflowdof
        zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ncount = ncount + 1
      zero_rows_local(ncount) = local_id*option%nflowdof
      zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
    enddo
#endif
  endif

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine createMphaseZeroArray


! ************************************************************************** !
!
! MphaseInitializeSolidReaction: Allocates and initializes arrays associated with
!                          mineral reactions
! author: Glenn Hammond
! date: 11/15/07
!
! ************************************************************************** !
subroutine MphaseInitializeSolidReaction(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscInt :: icell
  PetscReal, pointer :: phis_p(:)
  PetscErrorCode :: ierr
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  if (mphase_option%rk > 0.d0) then
    allocate(mphase_option%area_var(grid%nlmax))
    allocate(mphase_option%rate(grid%nlmax))
    call VecGetArrayF90(mphase_field%phis,phis_p,ierr)
    do icell = 1, grid%nlmax
      phis_p(icell) = mphase_option%phis0
      mphase_option%area_var(icell) = 1.d0
    enddo
    call VecRestoreArrayF90(mphase_field%phis,phis_p,ierr)
  endif
  
end subroutine MphaseInitializeSolidReaction

! ************************************************************************** !
!
! MphaseGetTecplotHeader: Returns a Tecplot file header
! author: 
! date: 
!
! ************************************************************************** !
function MphaseGetTecplotHeader(realization)

  use Realization_module
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: MphaseGetTecplotHeader
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = 'VARIABLES=' // &
           '"X [m]",' // &
           '"Y [m]",' // &
           '"Z [m]",' // &
           '"T [C]",' // &
           '"P [Pa]",' // &
           '"sl",' // &
           '"sg",' // &
           '"Ul",' // &
           '"Ug",'
  do i=1,option%nspec
    write(string2,'(''"Xl('',i2,'')",'')') i
    string = trim(string) // trim(string2)
  enddo
  do i=1,option%nspec
    write(string2,'(''"Xg('',i2,'')",'')') i
    string = trim(string) // trim(string2)
  enddo
  if (mphase_option%rk > 0.d0) then
    string = trim(string) // '"Volume Fraction"'
  endif
  string = trim(string) // ',"Phase"'
  if (associated(field%imat)) then
    string = trim(string) // ',"Material_ID"'
  endif
  
  MphaseGetTecplotHeader = string

end function MphaseGetTecplotHeader

! ************************************************************************** !
!
! MphaseGetVarFromArray: Extracts variables indexed by ivar and isubvar
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none
  
  PetscInt, parameter :: TEMPERATURE = 4
  PetscInt, parameter :: PRESSURE = 5
  PetscInt, parameter :: LIQUID_SATURATION = 6
  PetscInt, parameter :: GAS_SATURATION = 7
  PetscInt, parameter :: LIQUID_ENERGY = 8
  PetscInt, parameter :: GAS_ENERGY = 9
  PetscInt, parameter :: LIQUID_MOLE_FRACTION = 10
  PetscInt, parameter :: GAS_MOLE_FRACTION = 11
  PetscInt, parameter :: VOLUME_FRACTION = 12
  PetscInt, parameter :: PHASE = 13
  PetscInt, parameter :: MATERIAL_ID = 14

  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%grid
  field => realization%field

  call VecGetArrayF90(vec,vec_ptr,ierr)

  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)    
        select case(ivar)
          case(TEMPERATURE)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(PRESSURE)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(LIQUID_SATURATION)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(LIQUID_MOLE_FRACTION)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(LIQUID_ENERGY)
            vec_ptr(local_id) = 0.d0 ! to be provided
        end select
      enddo
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec2_ptr,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec2_ptr(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%iphas_loc,vec2_ptr,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = field%imat(grid%nL2G(local_id))
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine MphaseGetVarFromArray

end module MPHASE_module

