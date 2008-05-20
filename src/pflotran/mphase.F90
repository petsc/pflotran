module Mphase_module
  
  use Mphase_Aux_module
  
  implicit none
  
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
  PetscReal, parameter :: formeps   = 1.D-4
  PetscReal, parameter :: eps = 1.D-5 
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps   = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0  !1D-8
  

  PetscInt, parameter :; jh2o=1, jco2=2

  PetscReal, allocatable, save :: Resold_AR(:,:), Resold_FL(:,:)
  
  public MphaseResidual,MphaseJacobian, &
         MphaseUpdateFixedAccum,MphaseTimeCut,&
         MphaseSetup, MphaseNumericalJacTest, &
         MphaseGetVarFromArray, MphaseGetVarFromArrayAtCell, &
         MphaseMaxChange, MphaseUpdateSolution, &
         MphaseGetTecplotHeader, MphaseInitializeTimestep

contains

! ************************************************************************** !
!
! MphaseTimeCut: Resets arrays for time step cut
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscReal, pointer :: xx_p(:),yy_p(:)
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

  call VecCopy(field%flow_yy,field%flow_xx,ierr)
  call VecCopy(field%iphas_old,field%fiphas,ierr) 

end subroutine MphaseTimeCut

! ************************************************************************** !
!
! MphaseSetup: 
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseSetup(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MphaseSetup

! ************************************************************************** !
!
! MphaseSetupPatch: Creates arrays for auxilliary variables
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  type(Mphase_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Mphase => MphaseAuxCreate()
  
  ! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call MphaseAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%Mphase%aux_vars => aux_vars
  patch%aux%Mphase%num_aux = grid%ngmax
  
  allocate(delx(options%nflowdof, grid%ngmax))
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them  
  boundary_condition => patch%flow_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call MphaseAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%Mphase%aux_vars_bc => aux_vars_bc
  patch%aux%Mphase%num_aux_bc = sum_connection
  
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call MphaseCreateZeroArray(patch,option)

end subroutine MphaseSetupPatch

! ************************************************************************** !
! Mphaseinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
  function  MphaseInitGuessCheck(realization)
 
  use Realization_module
  use Level_module
  use Patch_module
  use Option_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: ipass, ipass0, ierr    

  option => realization%option
  cur_level => realization%level_list%first
  ipass = 1
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      ipass= MphaseInitGuessCheckPatch(realization)
      if(ipass<=0)then
        nullify(cur_level)
        exit 
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

   call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   if(option%commsize >1)then
      call MPI_ALLREDUCE(ipass,ipass0,ONE_INTEGER, MPI_INTEGER,MPI_SUM, &
           PETSC_COMM_WORLD,ierr)
      if(ipass0 < option%commsize) ipass=-1
   endif
   MphaseInitGuessCheck =ipass
 end function MphaseInitGuessCheck

! ************************************************************************** !
! Mphaseinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
 subroutine MPhaseUpdateReasonPatch(reason,realization)
   use Realization_module
   use Patch_module
   use Field_module
   use Option_module
   use Grid_module

  implicit none
 
  PetscInt, intent(out):: reason
  type(realization_type) :: realization  
  
  type(field_type), pointer :: field
  PetscReal, pointer :: xx_p(:),iphase_loc_p(:), yy_p(:) 
  PetscInt :: n,n0,re
  PetscInt :: re0, ierr, iipha
  
  option => realization%option
  field => realization%field  
  patch => realization%patch
  grid => patch%grid

  re=1
 
  if(re>0)then
     call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
     call VecGetArrayF90(field%yy, yy_p, ierr)
     call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); 
  
     do n = 1,grid%nlmax
!**** clu-Ignore inactive cells with inactive materials **************
        if (associated(patch%imat)) then
           if (patch%imat(grid%nL2G(n)) <= 0) cycle
        endif
        n0=(n-1)* option%nflowdof
        iipha=int(iphase_loc_p(grid%nL2G(n)))
  
! ******** Too huge change in pressure ****************     
        if(dabs(xx_p(n0 + 1)- yy_p(n0 + 1))> (10.0D0 * option%dpmxe))then
           re=0; print *,'huge change in p', xx_p(n0 + 1), yy_p(n0 + 1)
           exit
        endif

! ******** Too huge change in temperature ****************
        if(dabs(xx_p(n0 + 2)- yy_p(n0 + 2))> (10.0D0 * option%dtmpmxe))then
           re=0; print *,'huge change in T', xx_p(n0 + 2), yy_p(n0 + 2)
           exit
        endif
 
! ******* Check 0<=sat/con<=1 **************************
        select case(iipha)
        case (1)
           if(xx_p(n0 + 3) > 1.0D0)then
              re=0; exit
           endif
           if(xx_p(n0 + 3) < 0D0)then
              if(xx_p(n0 + 3) > -1D-3)then
                 xx_p(n0 + 3) =0.D0
              else  
                 re=0; exit
              endif
           endif
        case (2)
           if(xx_p(n0 + 3) > 1.0D0)then
              re=0; exit
           endif
           if(xx_p(n0 + 3) < 0D-0)then
              re=0; exit
           endif
        case (3)
           if(xx_p(n0 + 3) > 1.D0)then
              re=0; exit
           endif
           if(xx_p(n0 + 3) < 0.)then
              re=0; exit
           endif
        end select
     end do
  
     if(re<=0) print *,'Sat or Con out of Region at: ',n,iipha,xx_p(n0+1:n0+3)
     call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
     call VecRestoreArrayF90(grid%yy, yy_p, ierr)
     call VecRestoreArrayF90(grid%var, var_p, ierr) 
     call VecRestoreArrayF90(grid%iphas, iphase_p, ierr) 
  endif
  ! print *,' update reason', grid%myrank, re,n,grid%nlmax
 
  
end subroutine MPhaseUpdateReasonPatch


! ************************************************************************** !
!
! RichardsUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MPhaseUpdateReason(reason, update_realization)

  use Realization_module
  use Level_module
  use Patch_module
  implicit none

  type(realization_type) :: realization
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetScInt :: reason  

  reason = 1
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MPhaseUpdateReasonPatch(reason, realization)
        if(reason<=0)then
           nullify(cur_level)
           exit 
        endif
        cur_patch => cur_patch%next
     enddo
    cur_level => cur_level%next
 enddo

 call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  if(option%commsize >1)then
     call MPI_ALLREDUCE(re, re0,1, MPI_INTEGER,MPI_SUM, &
          PETSC_COMM_WORLD,ierr)
     if(re0<grid%commsize) re=0
  endif
  reason=re
  
  if(reason<=0 .and. option%myrank ==0) print *,'Sat or Con out of Region', re
end subroutine MPhaseUpdateReason

! ************************************************************************** !
! Mphaseinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
  function  MphaseInitGuessCheckPatch(realization)
   
     use span_wagner_module
     
    use Realization_module
    use Patch_module
    use Field_module
    use Grid_module
    use Option_module
    implicit none
    
    type(realization_type) :: realization

    type(grid_type), pointer :: grid
    type(patch_type), pointer :: patch
    type(option_type), pointer :: option
    type(field_type), pointer :: field
      
    PetscInt :: local_id
    PetscReal, pointer :: xx_p(:)


    patch => realization%patch
    grid => patch%grid
    option => realization%option
    field => realization%field
    
    call VecGetArrayF90(field%flow_xx,xx_p, ierr)
    
    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       !geh - Ignore inactive cells with inactive materials
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
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
    
    call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
    MphaseInitGuessCheckPatch = ierr
  end function MphaseInitGuessCheckPatch

! ************************************************************************** !
!
! MphaseUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Mphase problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseUpdateAuxVars(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MphaseUpdateAuxVars

! ************************************************************************** !
!
! MphaseUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                        the Mphase problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseUpdateAuxVarsPatch(realization)

  use Realization_module
  use Patch_module
  use Field_module
  use Option_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set
  type(Mphase_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  aux_vars => patch%aux%Mphase%aux_vars
  aux_vars_bc => patch%aux%Mphase%aux_vars_bc
  
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
   
    call MphaseAuxVarCompute(xx_loc_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_vars_elem(0), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       option)
    iphase_loc_p(ghosted_id) = iphase
  enddo

  boundary_condition => patch%flow_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      select case(boundary_condition%condition%itype(idof))
      case(DIRICHLET_BC)
         xxbc(:) = boundary_condition%aux_real_var(:,iconn)
      case(HYDROSTATIC_BC)
         xxbc(1) = boundary_condition%aux_real_var(1,iconn)
          xxbc(2:option%nflowdof) = &
               xx_loc_p((ghosted_id-1)*option%nflowdof+2:ghosted_id*option%nflowdof)
      case(CONST_TEMPERATURE)
   
      case(PRODUCTION_WELL) ! 102      
   
      case(NEUMANN_BC,ZERO_GRADIENT_BC)
         xxbc(:) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
      end select
      
      select case(boundary_condition%condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,SEEPAGE_BC)
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC, HYDROSTATIC_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      call MphaseAuxVarCompute(xxbc,aux_vars_bc(sum_connection)%aux_vars_elem(0), &
                         iphasebc, &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  patch%MphaseAux%aux_vars_up_to_date = .true.

end subroutine MphaseUpdateAuxVarsPatch

! ************************************************************************** !
!
! MphaseInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine MphaseInitializeTimestep(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call MphaseUpdateFixedAccumulation(realization)

end subroutine MphaseInitializeTimestep

! ************************************************************************** !
!
! MphaseUpdateSolution: Updates data in module after a successful time step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine MphaseUpdateSolution(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%flow_xx,realization%field%flow_yy,ierr)   
  call VecCopy(realization%field%iphas,realization%field%iphas_old,ierr)

! make room for hysteric s-Pc-kr

end subroutine MphaseUpdateSolution


! ************************************************************************** !
!
! MphaseUpdateFixedAccumulation: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !
subroutine MphaseUpdateFixedAccumulation(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MphaseUpdateFixedAccumulation

! ************************************************************************** !
!
! MphaseUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !
subroutine MphaseUpdateFixedAccumPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: aux_vars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  aux_vars => patch%aux%Mphase%aux_vars
    
  call VecGetArrayF90(field%flow_xx,xx_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
    call MphaseAuxVarCompute(xx_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_vars_elem(0), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       option)
    iphase_loc_p(ghosted_id) = iphase
    call MphaseAccumulation(aux_vars(ghosted_id)%aux_vars_elem(0), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)

#if 0
!  call MphaseNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine MphaseUpdateFixedAccumPatch

! ************************************************************************** !
!
! MphaseAccumulationDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Chuan Lu
! date: 12/13/07
!
! Only for analytical jacobian, not working yet
! ************************************************************************** !
#if 0
subroutine MphaseAccumulationDerivative(aux_var,por,vol,rock_dencpr,option, &
                                          sat_func,J)

  use Option_module
  use Material_module
  
  implicit none

  type(richards_auxvar_type) :: aux_var
  type(option_type) :: option
  PetscReal vol,por,rock_dencpr
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec !, iireac=1
  PetscReal :: porXvol, mol(option%nspec), eng

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: aux_var_pert
  PetscReal :: x(3), x_pert(3), pert, res(3), res_pert(3), J_pert(3,3)

  porXvol = por*vol
      
  J(1,1) = (aux_var%sat*aux_var%dden_dp+aux_var%dsat_dp*aux_var%den)*porXvol*aux_var%xmol(1)
  J(1,2) = (aux_var%sat*aux_var%dden_dt)*porXvol*aux_var%xmol(1)
  J(1,3) = 0.d0
  J(2,1) = (aux_var%sat*aux_var%dden_dp+aux_var%dsat_dp*aux_var%den)*porXvol*aux_var%xmol(2)
  J(2,2) = (aux_var%sat*aux_var%dden_dt)*porXvol*aux_var%xmol(2)
  J(2,3) = (aux_var%sat*aux_var%den)*porXvol
  J(3,1) = (aux_var%dsat_dp*aux_var%den*aux_var%u+ &
            aux_var%sat*aux_var%dden_dp*aux_var%u+ &
            aux_var%sat*aux_var%den*aux_var%du_dp)* &
           porXvol
  J(3,2) = aux_var%sat* &
           (aux_var%dden_dt*aux_var%u+ &  ! pull %sat outside
            aux_var%den*aux_var%du_dt)* &
           porXvol + (1.d0 - por)* vol * rock_dencpr 
  J(3,3) = 0.d0 

  if (option%numerical_derivatives) then
    allocate(aux_var_pert%xmol(option%nspec),aux_var_pert%diff(option%nspec))
    call MphaseAuxVarCopy(aux_var,aux_var_pert,option)
    x(1) = aux_var%pres
    x(2) = aux_var%temp
    x(3) = aux_var%xmol(2)
    call MphaseAccumulation(aux_var,por,vol,rock_dencpr,option,res)
    do ideriv = 1,3
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      x_pert(ideriv) = x_pert(ideriv) + pert
      call MphaseAuxVarCompute(x_pert,aux_var_pert,iphase,sat_func,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', aux_var%dvis_dp, (aux_var_pert%vis-aux_var%vis)/pert(ideriv)
!         print *, 'dkr_dp:', aux_var%dkr_dp, (aux_var_pert%kr-aux_var%kr)/pert(ideriv)
          print *, 'dsat_dp:', aux_var%dsat_dp, (aux_var_pert%sat-aux_var%sat)/pert
          print *, 'dden_dp:', aux_var%dden_dp, (aux_var_pert%den-aux_var%den)/pert
          print *, 'dkvr_dp:', aux_var%dkvr_dp, (aux_var_pert%kvr-aux_var%kvr)/pert
          print *, 'dh_dp:', aux_var%dh_dp, (aux_var_pert%h-aux_var%h)/pert
          print *, 'du_dp:', aux_var%du_dp, (aux_var_pert%u-aux_var%u)/pert
        case(2)
          print *, 'dden_dt:', aux_var%dden_dt, (aux_var_pert%den-aux_var%den)/pert
          print *, 'dkvr_dt:', aux_var%dkvr_dt, (aux_var_pert%kvr-aux_var%kvr)/pert
          print *, 'dh_dt:', aux_var%dh_dt, (aux_var_pert%h-aux_var%h)/pert
          print *, 'du_dt:', aux_var%du_dt, (aux_var_pert%u-aux_var%u)/pert
      end select     
#endif     
      call MphaseAccumulation(aux_var_pert,por,vol,rock_dencpr,option,res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo
    deallocate(aux_var_pert%xmol,aux_var_pert%diff)
    J = J_pert
  endif
   
end subroutine MphaseAccumulationDerivative

! ************************************************************************** !
!
! MphaseFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine MphaseFluxDerivative(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,sat_func_up,sat_func_dn,Jup,Jdn)
  use Option_module 
  use Material_module                             
  
  implicit none
  
  type(richards_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: v_darcy,area
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec
  PetscReal :: fluxm(option%nspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: ddifff_dp_up, ddifff_dp_dn, ddifff_dt_up, ddifff_dt_dn
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn, dden_ave_dt_up, dden_ave_dt_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn, dphi_dt_up, dphi_dt_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn, dukvr_dt_up, dukvr_dt_dn
  PetscReal :: duh_dp_up, duh_dp_dn, duh_dt_up, duh_dt_dn
  PetscReal :: dq_dp_up, dq_dp_dn, dq_dt_up, dq_dt_dn
  PetscReal :: duxmol_dxmol_up, duxmol_dxmol_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: aux_var_pert_up, aux_var_pert_dn
  PetscReal :: x_up(3), x_dn(3), x_pert_up(3), x_pert_dn(3), pert_up, pert_dn, &
            res(3), res_pert_up(3), res_pert_dn(3), J_pert_up(3,3), J_pert_dn(3,3)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0 
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dt_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddifff_dp_up = 0.d0
  ddifff_dp_dn = 0.d0
  ddifff_dt_up = 0.d0
  ddifff_dt_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_up = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_up = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_up = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_up = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_up = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_up = 0.d0
  duxmol_dxmol_dn = 0.d0
  
! Flow term
  if (aux_var_up%sat > sir_up .or. aux_var_dn%sat > sir_dn) then
    if (aux_var_up%sat <eps) then 
      upweight=0.d0
    else if (aux_var_dn%sat <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den
    dden_ave_dp_up = upweight*aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp
    dden_ave_dt_up = upweight*aux_var_up%dden_dt
    dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

    gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
              (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*aux_var_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

    dphi = aux_var_up%pres - aux_var_dn%pres  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
    dphi_dt_up = dgravity_dden_up*aux_var_up%dden_dt
    dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = aux_var_up%kvr
      dukvr_dp_up = aux_var_up%dkvr_dp
!      dukvr_dp_dn = 0.d0
      dukvr_dt_up = aux_var_up%dkvr_dt
!      dukvr_dt_dn = 0.d0
      
      uh = aux_var_up%h
      duh_dp_up = aux_var_up%dh_dp
!      duh_dp_dn = 0.d0
      duh_dt_up = aux_var_up%dh_dt
!      duh_dt_dn = 0.d0
      
      uxmol(1:option%nspec) = aux_var_up%xmol(1:option%nspec)
      duxmol_dxmol_up = 1.d0
!      duxmol_dxmol_dn = 0.d0
    else
      ukvr = aux_var_dn%kvr
!      dukvr_dp_up = 0.d0
      dukvr_dp_dn = aux_var_dn%dkvr_dp
!      dukvr_dt_up = 0.d0
      dukvr_dt_dn = aux_var_dn%dkvr_dt
      
      uh = aux_var_dn%h
!      duh_dp_up = 0.d0
      duh_dp_dn = aux_var_dn%dh_dp
!      duh_dt_up = 0.d0
      duh_dt_dn = aux_var_dn%dh_dt
      
      uxmol(1:option%nspec) = aux_var_dn%xmol(1:option%nspec)
!      duxmol_dxmol_up = 0.d0
      duxmol_dxmol_dn = 1.d0
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      dq_dt_up = Dq*(dukvr_dt_up*dphi+ukvr*dphi_dt_up)*area
      dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uxmol(1)
      Jup(1,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uxmol(1)
!      Jup(1,3:option%nflowdof) = 0d.0

      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(1)
      Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(1)
!      Jdn(1,3:option%nflowdof) = 0.d0
      do ispec=2,option%nspec 
        ! based on flux = q*density_ave*uxmol
        Jup(ispec,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uxmol(ispec)
        Jup(ispec,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uxmol(ispec)
        Jup(ispec,ispec+1) = q*density_ave*duxmol_dxmol_up

        Jdn(ispec,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(ispec)
        Jdn(ispec,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(ispec)
        Jdn(ispec,ispec+1) = q*density_ave*duxmol_dxmol_dn
      enddo
      ! based on flux = q*density_ave*uh
      Jup(option%nflowdof,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uh+q*density_ave*duh_dp_up
      Jup(option%nflowdof,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uh+q*density_ave*duh_dt_up
!      Jup(option%nflowdof,3:option%nflowdof) = 0d.0

      Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
      Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn
!      Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

    endif
  endif 
! Diffusion term   
! Note : average rule may not be correct  
  if ((aux_var_up%sat > eps) .and. (aux_var_dn%sat > eps)) then
!    difff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)* &
!                            (aux_var_up%den+aux_var_dn%den)
!    do ispec=1, option%nspec
!      fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
!                 (aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))* &
!                 (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
!    enddo 
    difff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)* &
                            (aux_var_up%den+aux_var_dn%den)
    ddifff_dp_up = diffdp * 0.25D0*(aux_var_up%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
                                    (aux_var_up%sat+aux_var_dn%sat)*aux_var_up%dden_dp)
    ddifff_dt_up = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_up%dden_dt
    ddifff_dp_dn = diffdp * 0.25D0*(aux_var_dn%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
                                    (aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dp)
    ddifff_dt_dn = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dt
                                    
    Jup(1,1) = Jup(1,1)+ddifff_dp_up*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    Jup(1,2) = Jup(1,2)+ddifff_dt_up*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
                                           
    Jdn(1,1) = Jdn(1,1)+ddifff_dp_dn*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    Jdn(1,2) = Jdn(1,2)+ddifff_dt_dn*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    do ispec=2, option%nspec
      Jup(ispec,1) = Jup(ispec,1)+ddifff_dp_up*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                                     (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jup(ispec,2) = Jup(ispec,2)+ddifff_dt_up*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                                     (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jup(ispec,ispec+1) = Jup(ispec,ispec+1)+difff*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))

      Jdn(ispec,1) = Jdn(ispec,1)+ddifff_dp_dn*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                             (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jdn(ispec,2) = Jdn(ispec,2)+ddifff_dt_dn*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                             (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jdn(ispec,ispec+1) = Jdn(ispec,ispec+1)+difff*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*(-1.d0)
    enddo  
  endif 

! conduction term
        
  Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
!  cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
  Jup(option%nflowdof,2) = Jup(option%nflowdof,2)+Dk*area
  Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2)+Dk*area*(-1.d0)
  Jup = Jup*option%flow_dt
  Jdn = Jdn*option%flow_dt
 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives) then
    allocate(aux_var_pert_up%xmol(option%nspec),aux_var_pert_up%diff(option%nspec))
    allocate(aux_var_pert_dn%xmol(option%nspec),aux_var_pert_dn%diff(option%nspec))
    call MphaseAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call MphaseAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_up(2) = aux_var_up%temp
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = aux_var_dn%pres
    x_dn(2) = aux_var_dn%temp
    x_dn(3) = aux_var_dn%xmol(2)
    call MphaseFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                      aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                      area,dist_gravity,upweight, &
                      option,v_darcy,res)
    do ideriv = 1,3
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn
      x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      call MphaseAuxVarCompute(x_pert_up,aux_var_pert_up,iphase,sat_func_up,option)
      call MphaseAuxVarCompute(x_pert_dn,aux_var_pert_dn,iphase,sat_func_dn,option)
      call MphaseFlux(aux_var_pert_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,res_pert_up)
      call MphaseFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_pert_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,res_pert_dn)
      J_pert_up(:,ideriv) = (res_pert_up(:)-res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_up%xmol,aux_var_pert_up%diff)
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jup = J_pert_up
    Jdn = J_pert_dn
  endif

end subroutine MphaseFluxDerivative

! ************************************************************************** !
!
! MphaseBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine MphaseBCFluxDerivative(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                                    area,dist_gravity,option, &
                                    sat_func_dn,Jdn)
  use Option_module
  use Material_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: area
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: v_darcy
  PetscReal :: fluxm(option%nspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: ddiff_dp_dn, ddiff_dt_dn
  PetscReal :: dden_ave_dp_dn, dden_ave_dt_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn, dphi_dt_dn
  PetscReal :: dukvr_dp_dn, dukvr_dt_dn
  PetscReal :: duh_dp_dn, duh_dt_dn
  PetscReal :: dq_dp_dn, dq_dt_dn
  PetscReal :: duxmol_dxmol_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: aux_var_pert_dn, aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(3), x_up(3), x_pert_dn(3), x_pert_up(3), pert_dn, res(3), &
            res_pert_dn(3), J_pert_dn(3,3)
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddiff_dp_dn = 0.d0
  ddiff_dt_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_dn = 0.d0
        
  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(RICHARDS_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (aux_var_up%sat > sir_dn .or. aux_var_dn%sat > sir_dn) then
        upweight=1.D0
        if (aux_var_up%sat < eps) then 
          upweight=0.d0
        else if (aux_var_dn%sat < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den
        dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp
        dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

        if (ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*aux_var_up%dden_dt
        endif
        
        gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
                  (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity
        
        dphi = aux_var_up%pres - aux_var_dn%pres + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
        dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt
        
        if (ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*aux_var_up%avgmw*dist_gravity*aux_var_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
          if (ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = aux_var_up%dkvr_dt
          endif
        else
          ukvr = aux_var_dn%kvr
          dukvr_dp_dn = aux_var_dn%dkvr_dp
          dukvr_dt_dn = aux_var_dn%dkvr_dt
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = aux_var_up%den
          if (ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dt_dn = aux_var_up%dden_dt
          endif
        else 
          density_ave = aux_var_dn%den
          dden_ave_dp_dn = aux_var_dn%dden_dp
          dden_ave_dt_dn = aux_var_dn%dden_dt
        endif 
        q = v_darcy * area
      endif

  end select

  if (v_darcy >= 0.D0) then
    uh = aux_var_up%h
    uxmol(:)=aux_var_up%xmol(1:option%nspec)
    if (ibndtype(RICHARDS_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = aux_var_up%dh_dp
    endif
    if (ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dt_dn = aux_var_up%dh_dt
    endif
    if (ibndtype(RICHARDS_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      duxmol_dxmol_dn = 1.d0
    endif
  else
    uh = aux_var_dn%h
    duh_dp_dn = aux_var_dn%dh_dp
    duh_dt_dn = aux_var_dn%dh_dt

    uxmol(:)=aux_var_dn%xmol(1:option%nspec)
    duxmol_dxmol_dn = 1.d0
  endif      

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(1)
  Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(1)
!  Jdn(1,3:option%nflowdof) = 0.d0
  do ispec=2,option%nspec 
    ! based on flux = q*density_ave*uxmol
    Jdn(ispec,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(ispec)
    Jdn(ispec,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(ispec)
    Jdn(ispec,ispec+1) = q*density_ave*duxmol_dxmol_dn
  enddo
      ! based on flux = q*density_ave*uh
  Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
  Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn
!  Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

  ! Diffusion term   
  select case(ibndtype(RICHARDS_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC) 
!      if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
!        diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
!        ddiff_dp_dn = diffdp * 0.25D0*(aux_var_dn%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
!                                      (aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dp)
!        ddiff_dt_dn = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dt
      if (aux_var_dn%sat > eps) then
        diff = diffdp * aux_var_dn%sat*aux_var_dn%den
        ddiff_dp_dn = diffdp * (aux_var_dn%dsat_dp*aux_var_dn%den+ &
                                aux_var_dn%sat*aux_var_dn%dden_dp)
        ddiff_dt_dn = diffdp * aux_var_dn%sat*aux_var_dn%dden_dt
                                    
        Jdn(1,1) = Jdn(1,1)+ddiff_dp_dn*aux_var_dn%diff(1)*&
                                        (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
        Jdn(1,2) = Jdn(1,2)+ddiff_dt_dn*aux_var_dn%diff(1)*&
                                        (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
        do ispec=2, option%nspec
          Jdn(ispec,1) = Jdn(ispec,1)+ddiff_dp_dn*aux_var_dn%diff(ispec)*&
                                                (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
          Jdn(ispec,2) = Jdn(ispec,2)+ddiff_dt_dn*aux_var_dn%diff(ispec)*&
                                                (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
          Jdn(ispec,ispec+1) = Jdn(ispec,ispec+1)+diff*aux_var_dn%diff(ispec)*(-1.d0)
        enddo  
      endif
  end select

  ! Conduction term
  select case(ibndtype(RICHARDS_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk =  Dk_dn / dd_up
      !cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
      Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2)+Dk*area*(-1.d0)
  end select

  Jdn = Jdn * option%flow_dt

  if (option%numerical_derivatives) then
    allocate(aux_var_pert_dn%xmol(option%nspec),aux_var_pert_dn%diff(option%nspec))
    allocate(aux_var_pert_up%xmol(option%nspec),aux_var_pert_up%diff(option%nspec))
    call MphaseAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call MphaseAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_up(2) = aux_var_up%temp
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = aux_var_dn%pres
    x_dn(2) = aux_var_dn%temp
    x_dn(3) = aux_var_dn%xmol(2)
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
    call MphaseBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                        por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                        area,dist_gravity,option,v_darcy,res)
    if (ibndtype(RICHARDS_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(RICHARDS_TEMPERATURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(RICHARDS_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    do ideriv = 1,3
      pert_dn = x_dn(ideriv)*perturbation_tolerance    
      x_pert_dn = x_dn
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif   
      call MphaseAuxVarCompute(x_pert_dn,aux_var_pert_dn,iphase,sat_func_dn,option)
      call MphaseAuxVarCompute(x_pert_up,aux_var_pert_up,iphase,sat_func_dn,option)
      call MphaseBCFlux(ibndtype,aux_vars,aux_var_pert_up,aux_var_pert_dn, &
                          por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                          area,dist_gravity,option,v_darcy,res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jdn = J_pert_dn
  endif

end subroutine MphaseBCFluxDerivative

#endif

! ************************************************************************** !
!
! MphaseAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !  
subroutine MphaseAccumulation(aux_var,por,vol,rock_dencpr,option,iireac,Res)

  use Option_module
  
  implicit none

  type(richards_auxvar_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
     
  PetscInt :: ispec, np, iireac
  PetscReal :: porXvol, mol(option%nspec), eng
  
 ! if (present(ireac)) iireac=ireac

  porXvol = por*vol
      
  mol=0.d0
  do np = 1, option%nphase
     do ispec=1, option%nspec  
        mol(ispec) = mol(ispec) + aux_var%sat(np) * &
             aux_var%den(np) * &
             aux_var%xmol(ispec + (np-1)*option%nspec)
     enddo
     eng = aux_var%sat(np) * 
     aux_var%den(np) * aux_var%u(np)
  enddo
  mol = mol * porXvol
  if(option%use_isoth == PETSC_FALSE) &
  eng = eng * porXvol + (1.d0 - por)* vol * rock_dencpr * aux_var%temp 
 
! Reaction terms here
! Note if iireac >0, then it is the node global index
  if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!H2O
     mol(1)= mol(1) - option%flow_dt * option%rtot(iireac,1)
     mol(2)= mol(2) - option%flow_dt * option%rtot(iireac,2)
  endif
  
   if(option%use_isoth)then
      Res(1:option%nflowdof)=mol(:)
   else
      Res(1:option%nflowdof-1)=mol(:)
      Res(option%nflowdof)=eng
   endif
end subroutine MphaseAccumulation

! ************************************************************************** !
!
! MphaseFlux: Computes the internal flux terms for the residual
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** ! 
subroutine MphaseFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(richards_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn(:)
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec
  PetscInt :: ispec
  PetscReal :: fluxm(option%nspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0  
  
! Flow term
  do np = 1, option%nphase
     if (aux_var_up%sat(np) > sir_up(np) .or. aux_var_dn%sat(np) > sir_dn(np)) then
        if (aux_var_up%sat(np) <eps) then 
           upweight=0.d0
        else if (aux_var_dn%sat(np) <eps) then 
           upweight=1.d0
        endif
        density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np) 
        
        gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
             (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
             * dist_gravity

        dphi = aux_var_up%pres - aux_var_dn%pres &
             - aux_var_up%pc(np) + aux_var_dn%pc(np) &
             + gravity

        ukvr=0.D0
        uh=0.D0
        uxmol=0.D0

        ! note uxmol only contains one phase xmol
        if (dphi>=0.D0) then
           ukvr = aux_var_up%kvr(np)
            if(option%use_isoth == PETSC_FALSE) uh = aux_var_up%h(np)
           uxmol(1:option%nspec) = aux_var_up%xmol((np-1)*option%nspec + 1 : np*option%nspec)
        else
           ukvr = aux_var_dn%kvr(np)
            if(option%use_isoth == PETSC_FALSE) uh = aux_var_dn%h(np)
           uxmol(1:option%nspec) = aux_var_dn%xmol((np-1)*option%nspec + 1 : np*option%nspec)
        endif
   

        if (ukvr>floweps) then
           v_darcy= Dq * ukvr * dphi
           vv_darcy(np)=v_darcy
           q = v_darcy * area
        
           do ispec=1, option%nspec 
              fluxm(ispec)=fluxm(ispec) + q * density_ave(np) * uxmol(ispec)
           enddo
           if(option%use_isoth == PETSC_FALSE) fluxe = fluxe + q*density_ave*uh 
        endif
     endif

! Diffusion term   
! Note : average rule may not be correct  
     if ((aux_var_up%sat(np) > eps) .and. (aux_var_dn%sat(np) > eps)) then
        difff = diffdp * 0.25D0*(aux_var_up%sat(np) + aux_var_dn%sat(np))* &
             (aux_var_up%den(np) + aux_var_dn%den(np))
        do ispec=1, option%nspec
           ind = ispec + (np-1)*option%nspec
           fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
                (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
        enddo
     endif
  enddo

! conduction term
  if(option%use_isoth == PETSC_FALSE) then     
     Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
     cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
     fluxe=fluxe + cond
  end if

  if(option%use_isoth)then
     Res(1:option%nflowdof) = fluxm(:) * option%flow_dt
  else
     Res(1:option%nflowdof-1) = fluxm(:) * option%flow_dt
     Res(option%nflowdof) = fluxe * option%flow_dt
  end if
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine MphaseFlux

! ************************************************************************** !
!
! MphaseBCFlux: Computes the  boundary flux terms for the residual
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !
subroutine MphaseBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(:)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm(option%nspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  do np = 1, option%nphase  
     select case(ibndtype(MPHASE_PRESSURE_DOF))
        ! figure out the direction of flow
     case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        Dq = perm_dn / dd_up
        ! Flow term
    
        if (aux_var_up%sat(np) > sir_dn(np) .or. aux_var_dn%sat(np) > sir_dn(np)) then
           upweight=1.D0
           if (aux_var_up%sat(np) < eps) then 
              upweight=0.d0
           else if (aux_var_dn%sat(np) < eps) then 
              upweight=1.d0
           endif
           density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np)
           
           gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
                (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
                * dist_gravity
       
           dphi = aux_var_up%pres - aux_var_dn%pres &
                - aux_var_up%pc(np) + aux_var_dn%pc(np) &
                + gravity
   
           if (dphi>=0.D0) then
              ukvr = aux_var_up%kvr(np)
           else
              ukvr = aux_var_dn%kvr(np)
           endif
     
           if (ukvr*Dq>floweps) then
              v_darcy = Dq * ukvr * dphi
           endif
        endif

     case(NEUMANN_BC)
        if (dabs(aux_vars(MPHASE_PRESSURE_DOF)) > floweps) then
           v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
           if (v_darcy > 0.d0) then 
              density_ave = aux_var_up%den(np)
           else 
              density_ave = aux_var_dn%den(np)
           endif
        endif

     end select

     q = v_darcy * area
     vv_darcy(np) = v_darcy
     uh=0.D0
     uxmol=0.D0
     
     if (v_darcy >= 0.D0) then
        if(option%use_isoth == PETSC_FALSE) uh = aux_var_up%h(np)
        uxmol(:)=aux_var_up%xmol((np-1)*option%nspec+1 : np * option%nspec)
     else
         if(option%use_isoth == PETSC_FALSE) uh = aux_var_dn%h(np)
        uxmol(:)=aux_var_dn%xmol((np-1)*option%nspec+1 : np * option%nspec)
     endif
    
     do ispec=1, option%nspec 
        fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
     enddo
      if(option%use_isoth == PETSC_FALSE) fluxe = fluxe + q*density_ave*uh
  enddo
     ! Diffusion term   
  select case(ibndtype(MPHASE_CONCENTRATION_DOF))
  case(DIRICHLET_BC) 
     !      if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
     !        diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
     if (aux_var_dn%sat(np) > eps) then
        diff = diffdp * aux_var_dn%sat(np) * aux_var_dn%den(np)
        do np = 1, option%nphase
           do ispec = 1, option%nspec
              fluxm(ispec) = fluxm(ispec) + diff * aux_var_dn%diff((np-1)* option%nspec+ispec)* &
                   (aux_var_up%xmol((np-1)* option%nspec+ispec) &
                   -aux_var_dn%xmol((np-1)* option%nspec+ispec))
           enddo
        enddo
     endif
  end select

  ! Conduction term
 if(option%use_isoth == PETSC_FALSE) then
    select case(ibndtype(MPHAASE_TEMPERATURE_DOF))
    case(DIRICHLET_BC, THERMAL_BC)
       Dk =  Dk_dn / dd_up
       cond = Dk*area*(aux_var_up%temp - aux_var_dn%temp) 
       fluxe=fluxe + cond
    end select
 end if

  Res(1:option%nspec)=fluxm(:)* option%flow_dt
  Res(option%nflowdof)=fluxe * option%flow_dt

end subroutine MphaseBCFlux

! ************************************************************************** !
!
! MphaseResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseResidual(snes,xx,r,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Field_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: ichange  

  field => realization%field
  grid => realization%patch%grid


!  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  
 ! check initial guess -----------------------------------------------
  ierr = MphaseInitGuessCheck(realization)

  ierr0 = 0
  if(option%commsize >1)then
    call MPI_ALLREDUCE(ierr, ierr0,ONE_INTEGER, MPI_INTEGER,MPI_SUM, PETSC_COMM_WORLD,ierr)
    if(ierr0 < option%commsize) then
      ierr=-1      
    endif
  endif
  
  if(ierr<0)then
    ierr = PETSC_ERR_ARG_DOMAIN
    if (option%myrank==0) print *,'table out of range: ',ierr0
    return
  endif 
  ! end check ---------------------------------------------------------


  ! Variable switching-------------------------------------------------
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseVarSwitchPatch((xx, realization, 0, ichange)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
! end switching ------------------------------------------------------

  ! Communication -----------------------------------------
  ! These 3 must be called before MphaseUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseResidualPatch(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MphaseResidual


! ************************************************************************** !
!
! MphaseVarSwitchPatch: Computes the residual equation at patch level
! author: Chuan Lu
! date: 3/10/08
!
! ************************************************************************** !

subroutine MphaseVarSwitchPatch(xx, realization, icri, ichange)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  
  use water_eos_module
  use gas_eos_module  
  use co2eos_module
  use span_wagner_module

  implicit none
  
  type(realization_type) :: realization
  
  Vec, intent(in) :: xx
  PetscInt :: icri,ichange 

  PetscReal, pointer :: xx_p(:), yy_p(:),iphase_loc_p(:)
  PetscScalar, pointer :: den(:)
  PetscInt :: ipr
  PetscInt :: iipha 
  PetscErrorCode :: ierr
! PetscInt :: index,i
  PetscReal :: p2,p,tmp,t
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  PetscReal :: ug,xphi,henry,sat_pressure
  PetscReal :: xmol(realization%option%nphase*realization%option%nspec),satu(realization%option%nphase)
! PetscReal :: xla,co2_poyn
  PetscInt :: local_id, ghosted_id, dof_offset
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
    
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  
! mphase code need assemble 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  den => patch%aux%Mphase%aux_vars%aux_var_elem(0)%den 
   
  ichange = 0   
  do local_id = 1,grid%nlmax
     ghosted_id = grid%nL2G(local_id)
     if (associated(grid%imat)) then
      if (grid%imat(ghosted_id) <= 0) cycle
    endif
     ipr=0 
     dof_offset=(local_id-1)* option%nflowdof
     iipha=iphase_loc_p(ghosted_id)
     p = xx_p(dof_offset+1)
     t= xx_p(dof_offset+2)
     select case(iipha) 
     case(1) 
        xmol(2)= xx_p(dof_offset+3)
        xmol(1)=1.D0 - xmol(2)
        satu(1)=1.D0; satu(2)=0.D0
     case(2)  
        xmol(4)= xx_p(dof_offset+3)
        xmol(3)=1.D0 - xmol(4)
        satu(1)=eps; satu(2)=1.D0
     case(3) 
        satu(2)= xx_p(dof_offset+3) 
        satu(1)= 1.D0- satu(2)
        xmol(3)= yh2o_in_co2; xmol(4)=1.D0-xmol(3)
     end select

! Pure CO2 phase properties ------------------------------------------    
     p2=p
!    p2=p*xmol(4)
    if(p2>=5d4)then
if(option%co2eos == 'EOS_SPAN_WAGNER')
 select case(grid%itable)  
    case(0,1,2,4,5)
      if( grid%itable >=4) then
        call co2_sw_interp(p2*1.D-6,t,dg,dddt,dddp,fg,&
        dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,grid%itable)
      else
        call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
        dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,grid%itable)
      endif
      dg= dg / fmwco2
      fg= fg * 1.D6 
      hg= hg * fmwco2
   ! Span-Wagner EOS with Bi-Cubic Spline interpolation
    case(3) 
     call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
      dg= dg / fmwco2
      fg= fg * 1.D6 
      hg= hg * fmwco2
  end select     
elseif(option%co2eos == 'EOS_MRK')
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
      call CO2( t,p2, dg,fg, xphi, hg)
      !call visco2(t,dg,visg)
      dg = dg / fmwco2
      hg = hg * fmwco2 *grid%scale
 endif
    else      
      call ideal_gaseos_noderiv(p2,t,grid%scale,dg,hg,ug)
      fg = p2
    endif
   
    xphi = fg/p2
    call PSAT(t, sat_pressure, ierr)
    sat_pressure =sat_pressure /1D5
    call Henry_duan_sun(t, p2 *1D-5, henry,xphi,option%m_nacl,option%m_nacl,sat_pressure)
    
    henry= 1.D8 / fmwh2o / henry / xphi !note: henry = H/phi
  
    wat_sat_x = sat_pressure*1.D5/p 
    co2_sat_x = (1.D0-wat_sat_x)/(henry/p-wat_sat_x)*henry/p  ! xmol(4) = xmol(2)*henry/p
!     tmp = 1.D0-tmp ! approximate form

select case(icri)
 case(0)
    select case(iipha)     
    case(1)
  
      xmol(4)=xmol(2)*henry/p   
      if (xmol(4) > 1.05D0*co2_sat_x) then
!      if (xmol(4) > 1.001D0*co2_sat_x .and. iipha==1) then
!     if (xmol(4) > (1.d0+1.d-6)*tmp .and. iipha==1) then
        write(*,'('' Liq -> 2ph '',''rank='',i6,'' it='',i3,'' n='',i8,'' p='',1pe10.4, &
      & '' T='',1pe10.4,'' Xl='',1pe11.4,'' xmol4='',1pe11.4, &
      & '' 1-Ps/P='',1pe11.4)') &
        grid%myrank,grid%iphch,n,xx_p(n0+1:n0+3),xmol(4),co2_sat_x
        iphase_p(n) = 3
        
        !Rachford-Rice initial guess: 1=H2O, 2=CO2
        k1 = wat_sat_x !sat_pressure*1.D5/p
        k2 = henry/p
        z1 = xmol(1); z2 = xmol(2)
        xg = ((1.d0-k2)*z2+(1.d0-k1)*z1)/((1.d0-k2)*(1.d0-k1)*(z1+z2))
        vmco2 = 1.d0/dg
        vmh2o = 1.D0 /den(1)   ! fmwh2o/0.9d3
       !calculate initial guess for sg
        sg = vmco2*xg/(vmco2*xg+vmh2o*(1.d0-xg))
        write(*,'(''Rachford-Rice: '',1p10e12.4)') k1,k2,z1,z2,xg,sg,den(1)
        xx_p(n0+3) = sg   
        ichange = 1
      endif

  
    case(2)   

      if (xmol(3) > wat_sat_x*1.05)then
!     if (xmol(3) > (1.d0+1.d-6)*tmp .and. iipha==2)then
        write(*,'('' Gas -> 2ph '',''rank='',i6,'' it='',i3,'' n='',i8, &
      & '' p= '',1pe10.4,'' T= '',1pe10.4,'' Xg= '',1pe11.4,'' Ps/P='', 1pe11.4)') &
        grid%myrank,grid%iphch,n,xx_p(n0+1:n0+3),wat_sat_x
        iphase_p(n) = 3
        xx_p(n0+3)=1.D0-formeps
        ichange = 1
      endif
 

    case(3)

      tmp = wat_sat_x
      xmol(2)= (1.D0-tmp)/(Henry/p-tmp) ! solve: x1+x2=1, y1+y2=1, y1=k1*x1, y2=k2*x2
!     xmol(2)= p*(1.D0-tmp)/Henry ! approximate form
      xmol(1)= 1.D0-xmol(2)
      xmol(3)= xmol(1)*tmp
      xmol(4)= 1.D0-xmol(3)
            
      if(satu(2)>=1.D0)then
        write(*,'('' 2ph -> Gas '',''rank='',i6,'' it='',i3,'' n='',i8, &
      & '' p='',1pe10.4,'' T='',1pe10.4,'' sg='',1p3e11.4)') &
        grid%myrank,grid%iphch,n,xx_p(n0+1:n0+3),satu(1),satu(2)
        iphase_p(n) = 2
        xx_p(n0 + 3) = 1.D0-1D-8
        ichange =1  
      else if(satu(2) <= 0.D0)then
        write(*,'('' 2ph -> Liq '',''rank= '',i6,'' it='',i3,'' n='',i8,'' p='',1pe10.4, &
      & '' T='',1pe10.4,'' sg ='',1pe11.4,'' sl='',1pe11.4,'' sg='',1pe11.4)')  &
        grid%myrank,grid%iphch,n,xx_p(n0+1:n0+3),satu(1),satu(2)
        iphase_p(n) = 1 ! 2ph -> Liq
        ichange = 1
        tmp = xmol(2) * 0.99
        xx_p(n0 + 3)=tmp
      endif


    end select

   case(1)
    
     select case(iipha)     
       case(2)   
         if (xmol(3) > wat_sat_x)then
           write(*,'(''** Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
         endif
       case(1) 
         xmol(4)=xmol(2)*henry/p 
         if (xmol(4) >co2_sat_x ) then
           write(*,'(''** Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3),xmol(4), co2_sat_x
         endif
       case(3) 
         if(satu(2)>1.D0 .and. iipha==3)then
            write(*,'(''** 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
         endif
         if(satu(2)<= 0.D0 .and. iipha==3)then
           write(*,'(''** 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3),satu(1),satu(2)
         endif
      end select
    end select

  enddo


  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%var, var_p,ierr)


end subroutine MphaseVarSwitchPatch
! ************************************************************************** !
!
! MphaseResidualPatch: Computes the residual equation at patch level
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseResidualPatch(snes,xx,r,realization,ierr)

  use water_eos_module

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  aux_vars => patch%aux%Mphase%aux_vars
  aux_vars_bc => patch%aux%Mphase%aux_vars_bc

 ! call MphaseUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MphaseAux%aux_vars_up_to_date = .false. 

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
 
  call VecGetArrayF90(field%flow_yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'

! Multiphase flash calculation is more expansive, so calculate once per iterration
#if 1
  ! Pertubations for aux terms --------------------------------
  
  do ng = 1, grid%ngmax
     if (associated(patch%imat)) then
        if (patch%imat(ng) <= 0) cycle
     endif
    
     istart =  (ng-1) * option%nflowdof +1 ; iend = istart -1 + option%nflowdof
     iphase =int(iphase_loc_p(ng))
 
     call MphaseAuxVarComputeNinc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0),iphase,&
          saturation_function,option)
     if (option%numerical_derivatives) then
        delx(1,ng) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac * 1.D-3
        delx(2,ng) = xx_loc_p((ng-1)*option%nflowdof+2)*dfac
        select case (iiphase)
        case (1)
           if(xx_loc_p((ng-1)*option%nflowdof+3) < 5D-5)then
               delx(3,ng) =  dfac*xx_loc_p((ng-1)*option%nflowdof+3)
           else
               delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
           endif
           if( delx(3,ng) < 1D-8 .and.  delx(3,ng)>=0.D0) delx(3,ng) =1D-8
           if( delx(3,ng) >-1D-8 .and.  delx(3,ng)<0.D0) delx(3,ng) =-1D-8
        case(2)  
           if(xx_loc_p((ng-1)*option%nflowdof+3) <0.9995)then
               delx(3,ng) =  dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
           else
               delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
           endif
           if( delx(3,ng) < 1D-8 .and.  delx(3,ng)>=0.D0) delx(3,ng) =1D-8
           if( delx(3,ng) >-1D-8 .and.  delx(3,ng)<0.D0) delx(3,ng) =-1D-8
        case(3)
           if(xx_loc_p((ng-1)*option%nflowdof+3) <=0.9)then
              delx(3,ng) = dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
           else
              delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
           endif
           
           if( delx(3,ng) < 1D-12 .and.  delx(3,ng)>=0.D0) delx(3,ng) = 1D-12
           if( delx(3,ng) >-1D-12 .and.  delx(3,ng)<0.D0) delx(3,ng) =-1D-12
        
           if(( delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3))>1.D0)then
              delx(3,ng) = (1.D0-xx_loc_p((ng-1)*option%nflowdof+3))*1D-6
           endif
           if(( delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3))<0.D0)then
              delx(3,ng) = xx_loc_p((ng-1)*option%nflowdof+3)*1D-6
           endif
        end select
 
        call MphaseAuxVarComputeWinc(xx_loc_p(istart:iend),delx(:,ng),&
             aux_vars(ng)%aux_var_elem(1:option%nflowdof),iphase, saturation_function,option)
     endif
  enddo
#endif


  Resold_AR=0.D0; ResOld_FL=0.D0; r_p = 0.d0

#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(0),porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
    ResAR_old(local_id, :)= Res(1:option%nflowdof)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%flow_source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%condition%num_sub_conditions > RICHARDS_CONCENTRATION_DOF) then
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
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (enthalpy_flag) then
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos_noderiv(tsrc1,aux_vars(ghosted_id)%aux_var_elem(0)%pres, &
                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o) &
                                               - qsrc1 *option%flow_dt
        if (enthalpy_flag) &
             r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - qsrc1*enth_src_h2o*option%flow_dt
      endif  
    
      if (csrc1 > 0.d0) then ! injection
!        call printErrMsg(option,"concentration source not yet implemented in Mphase")
      if(option%co2eos == 'EOS_SPAN_WAGNER')then
         !  span-wagner
            rho = aux_vars(ghosted_id)%aux_var_elem(0)%den(option%jco2)*option%fmwco2  
          select case(option%itable)  
            case(0,1,2,4,5)
              if( option%itable >=4) then
              call co2_sw_interp(aux_vars(ghosted_id)%aux_var_elem(0)%pres*1.D-6,&
                  tsrc1,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              else
              call co2_span_wagner(aux_vars(ghosted_id)%aux_var_elem(0)%pres*1.D-6,&
                  tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              endif 
             case(3) 
              call sw_prop(tsrc1,aux_vars(ghosted_id)%aux_var_elem(0)%pres*1.D-6,rho, &
                     enth_src_co2, eng, fg)
          end select     

         !  units: rho [kg/m^3]; csrc1 [kmol/s]
            enth_src_co2 = enth_src_co2 * option%fmwco2
      else if(option%co2eos == 'EOS_SPAN_WAGNER')then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]
            call CO2(tsrc1,aux_vars(ghosted_id)%aux_var_elem(0)%pres,&
            rho,fg, xphi,enth_src_co2)
            enth_src_co2 = enth_src_co2*option%fmwco2*option%scale
      else
         print *,'pflow mphase ERROR: Need specify CO2 EOS'
         STOP    
      endif
     
           
            r_p((n-1)*option%nflowdof + grid%jco2) = r_p((n-1)*option%nflowdof + option%jco2) - csrc1*option%flow_dt
            if (enthalpy_flag) &
              r_p(n*option%nflowdof) = r_p(n*option%nflowdof) - csrc1 * enth_src_co2 *option%flow_dt
            Resold_AR(n,grid%jco2)= Resold_AR(n,grid%jco2) - csrc1*option%flow_dt
            if (enthalpy_flag) &
              Resold_AR(n,option%nflowdof)= Resold_AR(n,option%nflowdof) - csrc1 * enth_src_co2*option%flow_dt

      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%flow_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = option%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  
! Then need fill up increments for BCs
    delxbc=0.D0;
    do idof =1, option%nflowdof   
       select case(boundary_condition%condition%itype(idof))
       case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof)=0.D0
       case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          iphasebc = int(iphase_loc_p(ng))
          delxbc=delx(1:grid%ndof,ghosted_id)
       end select
    enddo

    select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
    case(DIRICHLET_BC,SEEPAGE_BC)
       iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
    case(NEUMANN_BC,ZERO_GRADIENT_BC,HYDROSTATIC_BC)
       iphasebc=int(iphase_loc_p(ghosted_id))                               
    end select
 
    call MphaseAuxVarComputeNinc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),iphase,&
         saturation_function,option)

    call MphaseBCFlux(boundary_condition%condition%itype, &
         boundary_condition%aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection))%aux_var_elem(0), &
         aux_vars(ghosted_id)%)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         option%sir(1,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,Res)
    patch%boundary_velocities(:,sum_connection) = v_darcy(:)

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
      ResOld_AR(local_id,1:option%nflowdof) = ResOld_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = option%ckwet(ithrm_up)
      D_dn = option%ckwet(ithrm_dn)

      call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight,option,v_darcy,Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
      Resold_FL(sum_connection,1:option%nflowdof)= Res(1:option%nflowdof)
 
     if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif

! adjust residual to R/dt
  select case (option%idt_switch) 
  case(1) 
     r_p(:) = r_p(:)/option%flow_dt
  case(-1)
     if(option%flow_dt>1.D0) r_p(:) = r_p(:)/option%flow_dt
  end select
  
  do local_id = 1, grid%nlmax
     if (associated(patch%imat)) then
        if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
     endif

     istart = 1 + (n-1)*option%nflowdof
     if(volume_p(istart)>1.D0) r_p (istart:istart+2)=r_p(istart:istart+2)/volume_p(local_id)
     if(r_p(istart) >1E20 .or. r_p(istart) <-1E20) print *, r_p (istart:istart+2)
  enddo

! print *,'finished rp vol scale'
  if(option%use_isoth==PETSC_TRUE)then
     do local_id = 1, grid%nlmax  ! For each local node do...
        ghosted_id = grid%nL2G(local_id)   ! corresponding ghost index
        if (associated(grid%imat)) then
           if (grid%imat(ng) <= 0) cycle
        endif
        p1 = 3 + (n-1)*grid%ndof
        r_p(p1) = 0.D0 ! xx_loc_p(2 + (ng-1)*grid%ndof) - yy_p(p1-1)
     enddo
  endif
  !call VecRestoreArrayF90(r, r_p, ierr)

  
  if (n_zero_rows > 0) then
    do n=1,n_zero_rows
     r_p(zero_rows_local(n)) = 0.D0
    enddo
  endif


  if (patch%MphaseAux%inactive_cells_exist) then
    do i=1,patch%MphaseAux%n_zero_rows
      r_p(patch%MphaseAux%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Rresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Rxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine MphaseResidualPatch

! ************************************************************************** !
!
! MphaseJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Level_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseJacobianPatch(snes,xx,A,B,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MphaseJacobian

! ************************************************************************** !
!
! MphaseJacobianPatch: Computes the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine MphaseJacobianPatch(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,iphas_up,iphas_dn,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  PetscReal :: D_up, D_dn  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt ::  natural_id_up,natural_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
            Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  PetscReal :: delxbc(1:realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  
  PetscViewer :: viewer
  Vec :: debug_vec
!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  aux_vars => patch%aux%Mphase%aux_vars
  aux_vars_bc => patch%aux%Mphase%aux_vars_bc
  
! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
  flag = SAME_NONZERO_PATTERN

#if 0
!  call MphaseNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

 ResInc = 0.D0
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
     ghosted_id = grid%nL2G(local_id)
     !geh - Ignore inactive cells with inactive materials
     if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif
     iend = local_id*option%nflowdof
     istart = iend-option%nflowdof+1
     icap = int(icap_loc_p(ghosted_id))
     
     do nvar =1, option%nflowdof
        call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(nvar), &
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             option%dencpr(int(ithrm_loc_p(ghosted_id))), &
             option, res) 
        ResInc( local_id,:,nvar) =  ResInc(local_id,:,nvar) + Res(:)
     enddo
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%flow_source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%condition%num_sub_conditions > RICHARDS_CONCENTRATION_DOF) then
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

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
!      if (enthalpy_flag) then
!        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
!      endif         

      if (qsrc1 > 0.d0) then ! injection
       do nvar =1, option%nflowdof
         call wateos_noderiv(tsrc1,aux_vars(ghosted_id)%aux_var_elem(nvar)%pres,dw_kg,dw_mol,&
              enth_src_h2o,option%scale,ierr)        
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
          ! base on r_p() = r_p() - qsrc1*enth_src_h2o*option%flow_dt
           ResInc(local_id,option%jh2o,nvar)=  ResInc(local_id,option%jh2o,nvar) - qsrc1*option%flow_dt
           ResInc(local_id,option%nflowdof,nvar)=  ResInc(local_id,option%nflowdof,nvar) - qsrc1*enth_src_h2o*option%flow_dt
        enddo
     endif
    
      if (csrc1 > 0.d0) then ! injection
            do nvar=1,option%nflowdof     
              rho = aux_vars(ghosted_id)%aux_var_elem(nvar)%den(option%jco2)*option%fmwco2 
#ifdef EOS_SPAN_WAGNER
         !    span-wagner
           select case(option%itable)
             case(0,1,2, 4,5)
               if( option%itable >=4) then
                call co2_sw_interp(aux_vars(ghosted_id)%aux_var_elem(nvar)%pres*1.D-6, &
                    tsrc1,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              else
               call co2_span_wagner(aux_vars(ghosted_id)%aux_var_elem(nvar)%pres*1.D-6,&
                tsrc1+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
               endif
             case(3)
                call sw_prop(tsrc1,aux_vars(ghosted_id)%aux_var_elem(nvar)%pres*1.D-6,rho, &
                enth_src_co2, eng, fg)
          end select     
          enth_src_co2 = enth_src_co2 * grid%fmwco2  
#endif
#if EOS_MRK
! MRK eos [modified version from Kerrick and Jacobs (1981) and Weir et al. (1996).]     
              call CO2(tsrc1,aux_vars(ghosted_id)%aux_var_elem(nvar)%pres,&
                rho,fg, xphi,enth_src_co2 )
              enth_src_co2 = enth_src_co2 * grid%fmwco2 *option%scale
#endif

         !    units: rho [kg/m^3]; csrc1 [kmol/s]

              ResInc(local_id,option%jco2,nvar)=  ResInc(local_id,option%jco2,nvar) - csrc1*option%flow_dt
              ResInc(local_id,option%nflowdof,nvar)=  ResInc(local_id,option%nflowdof,nvar)&
                   - csrc1*enth_src_co2*option%flow_dt
           enddo
        endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
! Boundary conditions
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%flow_boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = option%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))

! Then need fill up increments for BCs
    delxbc=0.D0;
    do idof =1, option%nflowdof   
       select case(boundary_condition%condition%itype(idof))
       case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof)=0.D0
       case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          iphasebc = int(iphase_loc_p(ng))
          delxbc=delx(1:grid%ndof,ghosted_id)
       end select
    enddo

    select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
    case(DIRICHLET_BC,SEEPAGE_BC)
       iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
    case(NEUMANN_BC,ZERO_GRADIENT_BC,HYDROSTATIC_BC)
       iphasebc=int(iphase_loc_p(ghosted_id))                               
    end select
 
    call MphaseAuxVarComputeNinc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),iphase,&
         saturation_function,option)
    call MphaseAuxVarComputeWinc(xxbc,delxbc,&
         aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),iphase, saturation_function,option)
    
    do nvar=1,grid%ndof
       call MPHASERes_FLBCCont(boundary_condition%condition%itype, &
            aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
            aux_vars(ghosted_id)%aux_var_elem(nvar), &
            porosity_loc_p(ghosted_id), &
            tor_loc_p(ghosted_id), &
            option%sir(1,icap_dn), &
            cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
            cur_connection_set%area(iconn), &
            distance_gravity,option, &
            vv_darcy, Res)
       ResInc(local_id,1:grid%ndof,nvar) = ResInc(local_id,1:grid%ndof,nvar) - Res(1:grid%ndof)
    enddo
 enddo
    boundary_condition => boundary_condition%next
 enddo
#endif
! Set matrix values related to single node terms: Accumulation, Source/Sink, BC
  do local_id = 1, grid%nlmax  ! For each local node do...
     ghosted_id = grid%nL2G(local_id)
     !geh - Ignore inactive cells with inactive materials
     if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif

     ra=0.D0
     max_dev=0.D0
     do neq=1, grid%ndof
        do nvar=1, grid%ndof
           ra(neq,nvar)=(ResInc(local_id,neq,nvar)-ResOld_AR(local_id,neq))/delx(nvar,ghosted_id)
           if(max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
        enddo
     enddo

   select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%flow_dt
      case(-1)
        if(option%flow_dt>1) ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:) /option%flow_dt
    end select

     Jup=ra(1:grid%ndof,1:grid%ndof)
    
     if(volume_p(local_id)>1.D0 ) Jup=Jup / volume_p(local_id)
   
     ! if(n==1) print *,  blkmat11, volume_p(n), ra
     call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0    
  ResInc = 0.D0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or. &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      natural_id_up = grid%nG2N(ghosted_id_up)
      natural_id_dn = grid%nG2N(ghosted_id_dn)
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
    
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      iphas_up = iphase_loc_p(ghosted_id_up)
      iphas_dn = iphase_loc_p(ghosted_id_dn)

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      D_up = option%ckwet(ithrm_up)
      D_dn = option%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
      
      do nvar = 1, option%nflowdof 
         call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(1,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(1,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                           option, vv_darcy, )
            ra(:,nvar)= (Res(:)-ResOld_FL(nc,:))/grid%delx(nvar,ghosted_id_up)

         call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(1,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(1,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                           option, vv_darcy, )
         ra(:,nvar+grid%ndof)= (Res(:)-ResOld_FL(nc,:))/grid%delx(nvar,ghosted_id_down)
    enddo

    select case(option%idt_switch)
    case(1)
       ra =ra / option%flow_dt
    case(-1)  
       if(grid%dt>1)  ra =ra / option%flow_dt
    end select

 
    if (local_id_up > 0) then
       if(volume_p(local_id_up)>1.D0)then
          Jup= ra(:,1:option%ndof)/volume_p(local_id_up)
          jdn= ra(:, 1 + option%ndof:2 * option%ndof)/volume_p(local_id_up)
       endif
       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
    endif
    if (local_id_dn > 0) then
       if(volume_p(local_id_dn)>1.D0)then
          Jup= -ra(:,1:option%ndof)/volume_p(local_id_dn)
          jdn= -ra(:, 1 + option%ndof:2 * option%ndof)/volume_p(local_id_dn)
       endif
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
    endif
 enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero,ierr) 
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    elseif (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  if (patch%MphaseAux%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%MphaseAux%n_zero_rows, &
                          patch%MphaseAux%zero_rows_local_ghosted,f_up,ierr) 
  endif
#endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Rjacobian.out',viewer,ierr)
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

end subroutine MphaseJacobianPatch



! ************************************************************************** !
!
! RichardsCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine MphaseCreateZeroArray(patch,option)

  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  PetscInt :: flag = 0
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr
    
  grid => patch%grid
  
  n_zero_rows = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
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

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
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

  patch%aux%Richards%zero_rows_local => zero_rows_local
  patch%aux%Richards%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%Richards%n_zero_rows = n_zero_rows

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     PETSC_COMM_WORLD,ierr)
  if (flag > 0) patch%aux%Richards%inactive_cells_exist = .true.

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine MphaseCreateZeroArray

! ************************************************************************** !
!
! MphaseMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/15/08
!
! ************************************************************************** !
subroutine MphaseMaxChange(realization)

  use Realization_module
  use Level_module
  use Field_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal :: dcmac, dsmax, max_c, max_S  


  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  option%dpmax=0.D0
  optiondtmax=0.D0 
  option%dcmax=0.D0
  option%dsmax=0.D0
  dcmax=0.D0
  dsmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmax,ierr)

  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MphaseMaxChangePatch(realization, max_c, max_s)
      if(dcmax <max_c)  dcmax =max_c
      if(dsmax <max_s)  dsmax =max_s
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(grid%commsize >1)then
    call MPI_ALLREDUCE(dcmax, max_c,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dsmax, max_s,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    dcmax= max_C
    comp = max_s
  endif 
  option%dcmax=dcmax
  option%dsmax=dsmax

end subroutine MphaseMaxChange




subroutine MphaseMaxChangePatch(realization,  max_c, max_s)

  use Realization_module
  
  use Grid_module
  use Field_module
 
  implicit none
  type(realization_type) :: realization

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_p(:), yy_p(:),  iphase_loc_p(:),  iphase_old_loc_p(:) 
  PetscReal :: max_s, max_c 


  patch => realization%patch
  grid => patch%grid
  field => realization%field

  max_c=0.D0
  max_s=0.D0
  
 

  call VecGetArrayF90(field%flow_xx,xx_p, ierr)   
  call VecGetArrayF90(field%flow_yy,yy_p, ierr) 
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p, ierr)
  call VecGetArrayF90(field%iphas_old_loc,iphase_old_loc_p, ierr)

  do local_id =1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    
    if(int(iphase_loc_p(ghosted_id)) == int(iphase_old_loc_p(ghosted_id)))then
       cmp=dabs(xx_p(n0+3)-yy_p(n0+3))
       if(int(iphase_p(n))==1 .or.int(iphase_p(n))==2)then
          if(max_c<cmp) max_c = cmp
       endif
       if(int(iphase_p(n))==3)then
         if(max_s<cmp) max_s=cmp
      endif
   end if
  end do

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)   
  call VecRestoreArrayF90(field%flow_yy,yy_p, ierr) 
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_old_loc,iphase_old_loc_p, ierr)


  
end subroutine MphaseMaxChangePatch









! ************************************************************************** !
!
! RichardsLiteGetTecplotHeader: Returns Richards contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
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
  
  string = ',' // &
           '"T [C]",' // &
           '"P [Pa]",' // &
           '"PHASE",' // &
           '"Sg",' // &
           '"den(l),"'//&
           '"den(g)",' // &
           '"u(l),"'//&
           '"u(g)"' // &
   
  do i=1,option%nflowspec
    write(string2,'('',"Xl('',i2,'')"'')') i
    string = trim(string) // trim(string2)
  enddo
  
  MphaseGetTecplotHeader = string

end function MphaseGetTecplotHeader

! ************************************************************************** !
!
! RichardsGetVarFromArray: Extracts variables indexed by ivar and isubvar
!                          from Richards type
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine MphaseGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Field_module

  implicit none

  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: aux_vars(:)
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  if (.not.patch%aux%Mphase%aux_vars_up_to_date) call MphaseUpdateAuxVars(realization)

  aux_vars => patch%aux%Mphase%aux_vars
  
  call VecGetArrayF90(vec,vec_ptr,ierr)

  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)    
        select case(ivar)
          case(TEMPERATURE)
            vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%temp
          case(PRESSURE)
            vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%pres
          case(LIQUID_SATURATION)
            vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(1)
          case(LIQUID_DENSITY)
            vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%den(1)
          case(GAS_SATURATION)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(2)
          case(GAS_MOLE_FRACTION)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(4)
          case(GAS_ENERGY)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) ! still need implementation
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%den(2)
          case(LIQUID_MOLE_FRACTION)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2)
          case(LIQUID_ENERGY)
             vec_ptr(local_id) = aux_vars(ghosted_id)aux_var_elem(0)%%u(1)
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
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine MphaseGetVarFromArray

! ************************************************************************** !
!
! RichardsGetVarFromArrayAtCell: Returns variablesindexed by ivar, isubvar,
!                                 local id from Richards type
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function MphaseGetVarFromArrayAtCell(realization,ivar,isubvar,ghosted_id)

  use Realization_module
  use Grid_module
  use Patch_module
  use Option_module
  use Field_module

  implicit none

  PetscReal :: MphaseGetVarFromArrayAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: ghosted_id

  PetscReal :: value
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(richards_auxvar_type), pointer :: aux_vars(:)  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  if (.not.patch%aux%Mphase%aux_vars_up_to_date) call MphaseUpdateAuxVars(realization)

  aux_vars => patch%aux%Richards%aux_vars

  select case(ivar)
  case(TEMPERATURE)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%temp
  case(PRESSURE)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%pres
  case(LIQUID_SATURATION)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(1)
  case(LIQUID_DENSITY)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%den(1)
  case(GAS_SATURATION)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(2)
  case(GAS_MOLE_FRACTION)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(4)
  case(GAS_ENERGY)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
  case(GAS_DENSITY) ! still need implementation
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%den(2)
  case(LIQUID_MOLE_FRACTION)
     vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2)
  case(LIQUID_ENERGY)
     vec_ptr(local_id) = aux_vars(ghosted_id)aux_var_elem(0)%%u(1)
  end select
   
  RichardsGetVarFromArrayAtCell = value
  
end function MphaseGetVarFromArrayAtCell


! ************************************************************************** !
!
! RichardsDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MphaseDestroy(patch)

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call MphaseAuxDestroy(patch%aux%Richards)

end subroutine MphaseDestroy

! ************************************************************************** !
!
! MphaseCheckpointWrite: Writes vecs to checkpoint file
! author: 
! date: 
!
! ************************************************************************** !
subroutine MphaseCheckpointWrite(discretization, viewer)

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
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
subroutine MphaseCheckpointRead(discretization,viewer)

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call VecLoadIntoVector(viewer, global_var, ierr)
  call VecDestroy(global_var,ierr)
  ! solid volume fraction
  if (mphase_option%rk > 0.d0) then
    call VecLoadIntoVector(viewer, mphase_field%phis, ierr)
  endif  
  
end subroutine MphaseCheckpointRead



end module Mphase_module
