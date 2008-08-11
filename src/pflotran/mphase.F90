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
#include "include/finclude/petscerror.h"

! Cutoff parameters
  PetscReal, parameter :: formeps   = 1.D-4
  PetscReal, parameter :: eps = 1.D-5 
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps   = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0  !1D-8
  

  PetscInt, parameter :: jh2o=1, jco2=2

  PetscReal, allocatable, save :: Resold_AR(:,:), Resold_FL(:,:), delx(:,:)
  
  public MphaseResidual,MphaseJacobian, &
         MphaseUpdateFixedAccumulation,MphaseTimeCut,&
         MphaseSetup,MphaseUpdateReason,&
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
  call VecCopy(field%iphas_old_loc,field%iphas_loc,ierr) 

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
  use span_wagner_module
  use co2_sw_module
  use span_wagner_spline_module 
   
  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  if (realization%option%co2eos == EOS_SPAN_WAGNER)then
    select case(realization%option%itable)
       case(0,1,2)
         call initialize_span_wagner(realization%option%itable,realization%option%myrank)
       case(4,5)
         call initialize_span_wagner(0,realization%option%myrank)
         call initialize_sw_interp(realization%option%itable, realization%option%myrank)
       case(3)
         call sw_spline_read
       case default
         print *, 'Wrong table option : STOP'
      stop
    end select
  endif
 
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
  print *,' mph setup get patch'
  patch%aux%Mphase => MphaseAuxCreate()
  print *,' mph setup get Aux'
  ! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))
  print *,' mph setup get Aux alloc', grid%ngmax
  do ghosted_id = 1, grid%ngmax
    call MphaseAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%Mphase%aux_vars => aux_vars
  patch%aux%Mphase%num_aux = grid%ngmax
  print *,' mph setup get Aux init'

  allocate(delx(option%nflowdof, grid%ngmax))
  allocate(Resold_AR(grid%nlmax,option%nflowdof))
  allocate(Resold_FL(ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list),option%nflowdof))
  print *,' mph setup allocate app array'
   ! count the number of boundary connections and allocate
  ! aux_var data structures for them  
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(aux_vars_bc(sum_connection))
  print *,' mph setup get AuxBc alloc', sum_connection
  do iconn = 1, sum_connection
    call MphaseAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%Mphase%aux_vars_bc => aux_vars_bc
  patch%aux%Mphase%num_aux_bc = sum_connection
  option%numerical_derivatives = .true.

  print *,' mph setup get AuxBc point'
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
  
  PetscInt ::  MphaseInitGuessCheck
  type(realization_type) :: realization
  type(option_type), pointer:: option
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
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option 
  PetscReal, pointer :: xx_p(:),iphase_loc_p(:), yy_p(:) 
  PetscInt :: n,n0,re
  PetscInt :: re0, ierr, iipha
  
  option => realization%option
  field => realization%field  
  patch => realization%patch
  grid => patch%grid

  re=1
 
  if(re>0)then
     call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
     call VecGetArrayF90(field%flow_yy, yy_p, ierr)
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
    call VecRestoreArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr); 

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
subroutine MPhaseUpdateReason(reason, realization)

  use Realization_module
  use Level_module
  use Patch_module
  implicit none

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: reason

  PetscInt :: re, re0, ierr

  re = 1
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MPhaseUpdateReasonPatch(re, realization)
        if(re<=0)then
           nullify(cur_level)
           exit 
        endif
        cur_patch => cur_patch%next
     enddo
    cur_level => cur_level%next
 enddo

 call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  if(realization%option%commsize >1)then
     call MPI_ALLREDUCE(re, re0,1, MPI_INTEGER,MPI_SUM, &
          PETSC_COMM_WORLD,ierr)
     if(re0<realization%option%commsize) re=0
  endif
  reason=re
  
  if(reason<=0 .and. realization%option%myrank ==0) print *,'Sat or Con out of Region', re
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
    
    PetscInt :: MphaseInitGuessCheckPatch 
    type(realization_type) :: realization
    type(grid_type), pointer :: grid
    type(patch_type), pointer :: patch
    type(option_type), pointer :: option
    type(field_type), pointer :: field
      
    PetscInt :: local_id, ghosted_id, ierr, ipass
    PetscReal, pointer :: xx_p(:)


    patch => realization%patch
    grid => patch%grid
    option => realization%option
    field => realization%field
    
    call VecGetArrayF90(field%flow_xx,xx_p, ierr)
    
    ipass=1
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
          ipass=-1; exit  
       endif
       if(xx_p((local_id-1)*option%nflowdof+2)< t0_tab -273.15D0 &
            .or. xx_p((local_id-1)*option%nflowdof+2)>ntab_t*dt_tab + t0_tab-273.15D0)then
          ipass=-1; exit
       endif
    enddo

    call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
    MphaseInitGuessCheckPatch = ipass
  end function MphaseInitGuessCheckPatch

! ***********************************
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
  type(connection_set_type), pointer :: cur_connection_set
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
    if(.not. associated(realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr))then
       print*, 'error!!! saturation function not allocated', ghosted_id,icap_loc_p(ghosted_id)
    endif
   
    call MphaseAuxVarCompute_NINC(xx_loc_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_var_elem(0), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       realization%fluid_properties,option)
    iphase_loc_p(ghosted_id) = iphase
  enddo
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
    do idof=1,option%nflowdof
      select case(boundary_condition%flow_condition%itype(idof))
      case(DIRICHLET_BC)
         xxbc(:) = boundary_condition%flow_aux_real_var(:,iconn)
      case(HYDROSTATIC_BC)
         xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
          xxbc(2:option%nflowdof) = &
               xx_loc_p((ghosted_id-1)*option%nflowdof+2:ghosted_id*option%nflowdof)
     !  case(CONST_TEMPERATURE)
      
    !  case(PRODUCTION_WELL) ! 102      
   
      case(NEUMANN_BC,ZERO_GRADIENT_BC)
         xxbc(:) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
      end select
      enddo
      select case(boundary_condition%flow_condition%itype(1))
        case(DIRICHLET_BC,SEEPAGE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC, HYDROSTATIC_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      call MphaseAuxVarCompute_NINC(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0), &
                         iphasebc, &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         realization%fluid_properties, option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  patch%aux%Mphase%aux_vars_up_to_date = .true.

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
  call VecCopy(realization%field%iphas_loc,realization%field%iphas_old_loc,ierr)

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
  type(mphase_auxvar_type), pointer :: aux_vars(:)

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
    call MphaseAuxVarCompute_Ninc(xx_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_var_elem(0), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       realization%fluid_properties,option)
    iphase_loc_p(ghosted_id) = iphase
    call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(0), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,0, accum_p(istart:iend)) 
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
! MphaseAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !  
subroutine MphaseAccumulation(aux_var,por,vol,rock_dencpr,option,iireac,Res)

  use Option_module
  
  implicit none

  type(mphase_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
     
  PetscInt :: ispec, np, iireac
  PetscReal :: porXvol, mol(option%nflowspec), eng
  
 ! if (present(ireac)) iireac=ireac

  porXvol = por*vol
      
  mol=0.d0; eng=0.D0
  do np = 1, option%nphase
     do ispec=1, option%nflowspec  
        mol(ispec) = mol(ispec) + aux_var%sat(np) * &
             aux_var%den(np) * &
             aux_var%xmol(ispec + (np-1)*option%nflowspec)
     enddo
     eng = eng + aux_var%sat(np) * aux_var%den(np) * aux_var%u(np)
  enddo
  mol = mol * porXvol
 ! if(option%use_isoth == PETSC_FALSE) &
  eng = eng * porXvol + (1.d0 - por)* vol * rock_dencpr * aux_var%temp 
 
! Reaction terms here
! Note if iireac >0, then it is the node global index
 ! if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!H2O
 !    mol(1)= mol(1) - option%flow_dt * option%rtot(iireac,1)
 !    mol(2)= mol(2) - option%flow_dt * option%rtot(iireac,2)
 ! endif
  
   !if(option%use_isoth)then
   !   Res(1:option%nflowdof)=mol(:)
   !else
      Res(1:option%nflowdof-1)=mol(:)
      Res(option%nflowdof)=eng
  ! endif
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
  
  type(mphase_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),fluxe,q, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy =0.D0 
  
! Flow term
  do np = 1, option%nphase
     if (aux_var_up%sat(np) > sir_up(np) .or. aux_var_dn%sat(np) > sir_dn(np)) then
        upweight= dd_dn/(dd_up+dd_dn)
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

        v_darcy = 0.D0
        ukvr=0.D0
        uh=0.D0
        uxmol=0.D0

        ! note uxmol only contains one phase xmol
        if (dphi>=0.D0) then
           ukvr = aux_var_up%kvr(np)
           ! if(option%use_isoth == PETSC_FALSE)&
           uh = aux_var_up%h(np)
           uxmol(1:option%nflowspec) = aux_var_up%xmol((np-1)*option%nflowspec + 1 : np*option%nflowspec)
        else
           ukvr = aux_var_dn%kvr(np)
           ! if(option%use_isoth == PETSC_FALSE)&
           uh = aux_var_dn%h(np)
           uxmol(1:option%nflowspec) = aux_var_dn%xmol((np-1)*option%nflowspec + 1 : np*option%nflowspec)
        endif
   

        if (ukvr>floweps) then
           v_darcy= Dq * ukvr * dphi
           vv_darcy(np)=v_darcy
           q = v_darcy * area
        
           do ispec=1, option%nflowspec 
              fluxm(ispec)=fluxm(ispec) + q * density_ave * uxmol(ispec)
           enddo
          ! if(option%use_isoth == PETSC_FALSE)&
            fluxe = fluxe + q*density_ave*uh 
        endif
     endif

! Diffusion term   
! Note : average rule may not be correct  
     if ((aux_var_up%sat(np) > eps) .and. (aux_var_dn%sat(np) > eps)) then
        difff = diffdp * 0.25D0*(aux_var_up%sat(np) + aux_var_dn%sat(np))* &
             (aux_var_up%den(np) + aux_var_dn%den(np))
        do ispec=1, option%nflowspec
           ind = ispec + (np-1)*option%nflowspec
           fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
                (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
        enddo
     endif
  enddo

! conduction term
  !if(option%use_isoth == PETSC_FALSE) then     
     Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
     cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
     fluxe=fluxe + cond
 ! end if

  !if(option%use_isoth)then
  !   Res(1:option%nflowdof) = fluxm(:) * option%flow_dt
 ! else
     Res(1:option%nflowdof-1) = fluxm(:) * option%flow_dt
     Res(option%nflowdof) = fluxe * option%flow_dt
 ! end if
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
  type(mphase_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(:)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  do np = 1, option%nphase  
     select case(ibndtype(1))
        ! figure out the direction of flow
     case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        Dq = perm_dn / dd_up
        ! Flow term
        ukvr=0.D0
        v_darcy=0.D0 
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
        v_darcy = 0.D0
        if (dabs(aux_vars(1)) > floweps) then
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
        !if(option%use_isoth == PETSC_FALSE)&
         uh = aux_var_up%h(np)
        uxmol(:)=aux_var_up%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     else
         !if(option%use_isoth == PETSC_FALSE)&
        uh = aux_var_dn%h(np)
        uxmol(:)=aux_var_dn%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     endif
    
     do ispec=1, option%nflowspec 
        fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
     enddo
      !if(option%use_isoth == PETSC_FALSE) &
      fluxe = fluxe + q*density_ave*uh
 !print *,'FLBC', ibndtype(1),np, ukvr, v_darcy, uh, uxmol
   enddo
     ! Diffusion term   
  select case(ibndtype(3))
  case(DIRICHLET_BC) 
     !      if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
     !        diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
        do np = 1, option%nphase
          if(aux_var_up%sat(np)>eps .and. aux_var_dn%sat(np)>eps)then
              diff =diffdp * 0.25D0*(aux_var_up%sat(np)+aux_var_dn%sat(np))*&
                    (aux_var_up%den(np)+aux_var_up%den(np))
           do ispec = 1, option%nflowspec
              fluxm(ispec) = fluxm(ispec) + diff * aux_var_dn%diff((np-1)* option%nflowspec+ispec)* &
                   (aux_var_up%xmol((np-1)* option%nflowspec+ispec) &
                   -aux_var_dn%xmol((np-1)* option%nflowspec+ispec))
           enddo
          endif         
        enddo
     
  end select

  ! Conduction term
! if(option%use_isoth == PETSC_FALSE) then
    select case(ibndtype(2))
    case(DIRICHLET_BC, 4)
       Dk =  Dk_dn / dd_up
       cond = Dk*area*(aux_var_up%temp - aux_var_dn%temp) 
       fluxe=fluxe + cond
    end select
! end if

  Res(1:option%nflowspec)=fluxm(:)* option%flow_dt
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
  use Discretization_module
  use Field_module
  use Option_module
  use grid_module 

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: ichange  

  field => realization%field
  grid => realization%patch%grid
  option => realization%option
  discretization => realization%discretization
  
 
!  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
 ! check initial guess -----------------------------------------------
  ierr = MphaseInitGuessCheck(realization)
  if(ierr<0)then
    !ierr = PETSC_ERR_ARG_OUTOFRANGE
    if (option%myrank==0) print *,'table out of range: ',ierr
    call SNESSetFunctionDomainError() 
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
      call MphaseVarSwitchPatch(xx, realization, 0, ichange)
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
  use span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use span_wagner_module

  implicit none
  
  type(realization_type) :: realization
  
  Vec, intent(in) :: xx
  PetscInt :: icri,ichange 

  PetscReal, pointer :: xx_p(:), yy_p(:),iphase_loc_p(:)
  PetscReal :: den(realization%option%nphase)
  PetscInt :: ipr
  PetscInt :: iipha 
  PetscErrorCode :: ierr
! PetscInt :: index,i
  PetscReal :: p2,p,tmp,t
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  PetscReal :: ug,xphi,henry,sat_pressure
  PetscReal :: k1, k2, z1, z2,xg, vmco2, vmh2o,sg
  PetscReal :: xmol(realization%option%nphase*realization%option%nflowspec),&
               satu(realization%option%nphase)
  PetscReal :: yh2o_in_co2, wat_sat_x, co2_sat_x
! PetscReal :: xla,co2_poyn
  PetscInt :: local_id, ghosted_id, dof_offset
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  
  patch => realization%patch  
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  
! mphase code need assemble 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  
  ichange = 0   
  do local_id = 1,grid%nlmax
     ghosted_id = grid%nL2G(local_id)
     if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
     ipr=0 
     dof_offset=(local_id-1)* option%nflowdof
     iipha=iphase_loc_p(ghosted_id)
     p = xx_p(dof_offset+1)
     t= xx_p(dof_offset+2)
    den(1:option%nphase) = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(1:option%nphase)
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
       if(option%co2eos == EOS_SPAN_WAGNER)then
         select case(option%itable)  
          case(0,1,2,4,5)
          if(option%itable >=4) then
            call co2_sw_interp(p2*1.D-6,t,dg,dddt,dddp,fg,&
                dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
           else
             call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
                dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
          endif
          dg= dg / option%fmwco2
          fg= fg * 1.D6 
          hg= hg * option%fmwco2
! Span-Wagner EOS with Bi-Cubic Spline interpolation
          case(3) 
            call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
            dg= dg / option%fmwco2
            fg= fg * 1.D6 
            hg= hg * option%fmwco2
          end select     
        elseif(option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
          call CO2( t,p2, dg,fg, xphi, hg)
          dg = dg / option%fmwco2
          hg = hg * option%fmwco2 *option%scale
       endif
    else      
       call ideal_gaseos_noderiv(p2,t,option%scale,dg,hg,ug)
       fg = p2
    endif
   
    xphi = fg/p2
    call PSAT(t, sat_pressure, ierr)
    sat_pressure =sat_pressure /1D5
    call Henry_duan_sun(t, p2 *1D-5, henry,xphi,option%m_nacl,option%m_nacl,sat_pressure)
    
    henry= 1.D8 / option%fmwh2o / henry / xphi !note: henry = H/phi
  
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
        write(*,'('' Liq -> 2ph '',''rank='',i6,'' n='',i8,'' p='',1pe10.4, &
      & '' T='',1pe10.4,'' Xl='',1pe11.4,'' xmol4='',1pe11.4, &
      & '' 1-Ps/P='',1pe11.4)') &
        option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3),xmol(4),co2_sat_x
        iphase_loc_p(ghosted_id) = 3
        
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
        xx_p(dof_offset+3) = sg   
        ichange = 1
      endif

  
    case(2)   

      if (xmol(3) > wat_sat_x*1.05)then
!     if (xmol(3) > (1.d0+1.d-6)*tmp .and. iipha==2)then
        write(*,'('' Gas -> 2ph '',''rank='',i6,'' n='',i8, &
      & '' p= '',1pe10.4,'' T= '',1pe10.4,'' Xg= '',1pe11.4,'' Ps/P='', 1pe11.4)') &
        option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3),wat_sat_x
        iphase_loc_p(ghosted_id) = 3
        xx_p(dof_offset+3)=1.D0-formeps
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
        write(*,'('' 2ph -> Gas '',''rank='',i6,'' n='',i8, &
      & '' p='',1pe10.4,'' T='',1pe10.4,'' sg='',1p3e11.4)') &
        option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
        iphase_loc_p(ghosted_id) = 2
        xx_p(dof_offset+3) = 1.D0-1D-8
        ichange =1  
      else if(satu(2) <= 0.D0)then
        write(*,'('' 2ph -> Liq '',''rank= '',i6,'' n='',i8,'' p='',1pe10.4, &
      & '' T='',1pe10.4,'' sg ='',1pe11.4,'' sl='',1pe11.4,'' sg='',1pe11.4)')  &
        option%myrank,local_id, xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
        iphase_loc_p(ghosted_id) = 1 ! 2ph -> Liq
        ichange = 1
        tmp = xmol(2) * 0.99
        xx_p(dof_offset+3)=tmp
      endif


    end select

   case(1)
    
     select case(iipha)     
       case(2)   
         if (xmol(3) > wat_sat_x)then
           write(*,'(''** Gas -> 2ph '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3)
         endif
       case(1) 
         xmol(4)=xmol(2)*henry/p 
         if (xmol(4) >co2_sat_x ) then
           write(*,'(''** Liq -> 2ph '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3),xmol(4), co2_sat_x
         endif
       case(3) 
         if(satu(2)>1.D0 .and. iipha==3)then
            write(*,'(''** 2ph -> Gas '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3)
         endif
         if(satu(2)<= 0.D0 .and. iipha==3)then
           write(*,'(''** 2ph -> Liq '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
         endif
      end select
    end select
   enddo


  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p,ierr)

end subroutine MphaseVarSwitchPatch
! ************************************************************************** !
!
! MphaseResidualPatch: Computes the residual equation at patch level
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseResidualPatch(snes,xx,r,realization,ierr)

  use Connection_module
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
   use water_eos_module
  use gas_eos_module  
  use co2eos_module
  use span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use span_wagner_module

  
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
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(mphase_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
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
 
 
! Multiphase flash calculation is more expansive, so calculate once per iterration
#if 1
  ! Pertubations for aux terms --------------------------------
  do ng = 1, grid%ngmax
     if(grid%nG2L(ng)<0)cycle
     if (associated(patch%imat)) then
        if (patch%imat(ng) <= 0) cycle
     endif
        
     istart =  (ng-1) * option%nflowdof +1 ; iend = istart -1 + option%nflowdof
     iphase =int(iphase_loc_p(ng))
     call MphaseAuxVarCompute_Ninc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0),iphase,&
          realization%saturation_function_array(int(icap_loc_p(ng)))%ptr,&
          realization%fluid_properties,option)

     if (option%numerical_derivatives) then
        delx(1,ng) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac * 1.D-3
        delx(2,ng) = xx_loc_p((ng-1)*option%nflowdof+2)*dfac
        select case (iphase)
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
       call MphaseAuxVarCompute_Winc(xx_loc_p(istart:iend),delx(:,ng),&
            aux_vars(ng)%aux_var_elem(1:option%nflowdof),iphase,&
            realization%saturation_function_array(int(icap_loc_p(ng)))%ptr,&
            realization%fluid_properties,option)
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
                              option,1,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
    !print *,'REs, acm: ', res
    Resold_AR(local_id, :)= Res(1:option%nflowdof)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    !print *, 'RES s/s begin'
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = .true.
   ! else
   !   enthalpy_flag = .false.
   ! endif
      

    qsrc1 = source_sink%flow_condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
    !if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)
    hsrc1=0D0
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
       ! print *,'REs: S/S:',iconn, local_id, qsrc1,tsrc1,csrc1 
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
      if(option%co2eos == EOS_SPAN_WAGNER)then
         !  span-wagner
            rho = aux_vars(ghosted_id)%aux_var_elem(0)%den(jco2)*option%fmwco2  
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
      else if(option%co2eos == EOS_SPAN_WAGNER)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]
            call CO2(tsrc1,aux_vars(ghosted_id)%aux_var_elem(0)%pres,&
            rho,fg, xphi,enth_src_co2)
            enth_src_co2 = enth_src_co2*option%fmwco2*option%scale
      else
         print *,'pflow mphase ERROR: Need specify CO2 EOS'
         STOP    
      endif
     
           
            r_p((local_id-1)*option%nflowdof + jco2) = r_p((local_id-1)*option%nflowdof + jco2) - csrc1*option%flow_dt
            if (enthalpy_flag) &
              r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - csrc1 * enth_src_co2 *option%flow_dt
            Resold_AR(local_id,jco2)= Resold_AR(local_id,jco2) - csrc1*option%flow_dt
            if (enthalpy_flag) &
              Resold_AR(local_id,option%nflowdof)= Resold_AR(local_id,option%nflowdof) - csrc1 * enth_src_co2*option%flow_dt

      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
        
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
    do idof =1, option%nflowdof   
       select case(boundary_condition%flow_condition%itype(idof))
       case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
       case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          iphase = int(iphase_loc_p(ghosted_id))
       end select
    enddo

    select case(boundary_condition%flow_condition%itype(3))
    case(DIRICHLET_BC,SEEPAGE_BC)
       iphase = boundary_condition%flow_aux_int_var(1,iconn)
    case(NEUMANN_BC,ZERO_GRADIENT_BC,HYDROSTATIC_BC)
       iphase=int(iphase_loc_p(ghosted_id))                               
    end select
 
    call MphaseAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),iphase,&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)

    call MphaseBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(0), &
         aux_vars(ghosted_id)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         option%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,Res)
    patch%boundary_velocities(:,sum_connection) = v_darcy(:)
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
    Resold_AR(local_id,1:option%nflowdof) = ResOld_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
  enddo
  boundary_condition => boundary_condition%next
 enddo
#endif
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
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

      patch%internal_velocities(:,sum_connection) = v_darcy(:)
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

     istart = 1 + (local_id-1)*option%nflowdof
     if(volume_p(local_id)>1.D0) r_p (istart:istart+2)=r_p(istart:istart+2)/volume_p(local_id)
     if(r_p(istart) >1E20 .or. r_p(istart) <-1E20) print *, r_p (istart:istart+2)
  enddo

! print *,'finished rp vol scale'
  if(option%use_isoth) then
     do local_id = 1, grid%nlmax  ! For each local node do...
        ghosted_id = grid%nL2G(local_id)   ! corresponding ghost index
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) <= 0) cycle
        endif
        istart = 3 + (local_id-1)*option%nflowdof
        r_p(istart) = 0.D0 ! xx_loc_p(2 + (ng-1)*option%nflowdof) - yy_p(p1-1)
     enddo
  endif
  !call VecRestoreArrayF90(r, r_p, ierr)


  if (patch%aux%Mphase%inactive_cells_exist) then
    do i=1,patch%aux%Mphase%n_zero_rows
      r_p(patch%aux%Mphase%zero_rows_local(i)) = 0.d0
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
! author: Chuan Lu
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
! author: Chuan Lu
! date: 12/13/07
!
! ************************************************************************** !
subroutine MphaseJacobianPatch(snes,xx,A,B,flag,realization,ierr)

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  use water_eos_module
  use gas_eos_module  
  use co2eos_module
  use span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use span_wagner_module
    
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
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof,&
           realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(mphase_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  
  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: iphasebc                
  
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
             option,1, res) 
        ResInc( local_id,:,nvar) =  ResInc(local_id,:,nvar) + Res(:)
     enddo
     
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = .true.
   ! else
   !   enthalpy_flag = .false.
   ! endif

    qsrc1 = source_sink%flow_condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection_set
    
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
           ResInc(local_id,jh2o,nvar)=  ResInc(local_id,jh2o,nvar) - qsrc1*option%flow_dt
          if (enthalpy_flag) & 
           ResInc(local_id,option%nflowdof,nvar)=  ResInc(local_id,option%nflowdof,nvar) - qsrc1*enth_src_h2o*option%flow_dt
        enddo
     endif
    
      if (csrc1 > 0.d0) then ! injection
            do nvar=1,option%nflowdof     
              rho = aux_vars(ghosted_id)%aux_var_elem(nvar)%den(jco2)*option%fmwco2 
            if(option%co2eos == EOS_SPAN_WAGNER)then
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
          enth_src_co2 = enth_src_co2 * option%fmwco2  
          elseif(option%co2eos == EOS_MRK)then 
! MRK eos [modified version from Kerrick and Jacobs (1981) and Weir et al. (1996).]     
              call CO2(tsrc1,aux_vars(ghosted_id)%aux_var_elem(nvar)%pres,&
                rho,fg, xphi,enth_src_co2 )
              enth_src_co2 = enth_src_co2 * option%fmwco2 *option%scale
          endif

         !    units: rho [kg/m^3]; csrc1 [kmol/s]

              ResInc(local_id,jco2,nvar)=  ResInc(local_id,jco2,nvar) - csrc1*option%flow_dt
              if (enthalpy_flag) &
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
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
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
       select case(boundary_condition%flow_condition%itype(idof))
       case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof)=0.D0
       case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          !iphasebc = int(iphase_loc_p(ghosted_id))
          delxbc(idof)=delx(idof,ghosted_id)
       end select
    enddo
    !print *,'BC:',boundary_condition%flow_condition%itype, xxbc, delxbc


    select case(boundary_condition%flow_condition%itype(3))
    case(DIRICHLET_BC,SEEPAGE_BC)
       iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
    case(NEUMANN_BC,ZERO_GRADIENT_BC,HYDROSTATIC_BC)
       iphasebc=int(iphase_loc_p(ghosted_id))                               
    end select
 
    call MphaseAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),iphasebc,&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)
    call MphaseAuxVarCompute_Winc(xxbc,delxbc,&
         aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),iphasebc,&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties,option)
    
    do nvar=1,option%nflowdof
       call MphaseBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
         aux_vars(ghosted_id)%aux_var_elem(nvar), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         option%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         vv_darcy,Res)
       ResInc(local_id,1:option%nflowdof,nvar) = ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
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
     do neq=1, option%nflowdof
        do nvar=1, option%nflowdof
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

     Jup=ra(1:option%nflowdof,1:option%nflowdof)
     if(volume_p(local_id)>1.D0 ) Jup=Jup / volume_p(local_id)
   
     ! if(n==1) print *,  blkmat11, volume_p(n), ra
     call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
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
     ! natural_id_up = grid%nG2N(ghosted_id_up)
     ! natural_id_dn = grid%nG2N(ghosted_id_dn)
   
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
                          tor_loc_p(ghosted_id_up),option%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
            ra(:,nvar)= (Res(:)-ResOld_FL(iconn,:))/delx(nvar,ghosted_id_up)

         call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_dn),&
                          tor_loc_p(ghosted_id_dn),option%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
         ra(:,nvar+option%nflowdof)= (Res(:)-ResOld_FL(iconn,:))/delx(nvar,ghosted_id_dn)
    enddo

    select case(option%idt_switch)
    case(1)
       ra =ra / option%flow_dt
    case(-1)  
       if(option%flow_dt>1)  ra =ra / option%flow_dt
    end select
    
    if (local_id_up > 0) then
       voltemp=1.D0
       if(volume_p(local_id_up)>1.D0)then
         voltemp = 1.D0/volume_p(local_id_up)
       endif
       Jup(:,1:option%nflowdof)= ra(:,1:option%nflowdof)*voltemp !11
       jdn(:,1:option%nflowdof)= ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !12

       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
    endif
    if (local_id_dn > 0) then
       voltemp=1.D0
       if(volume_p(local_id_dn)>1.D0)then
         voltemp=1.D0/volume_p(local_id_dn)
       endif
       Jup(:,1:option%nflowdof)= -ra(:,1:option%nflowdof)*voltemp !21
       jdn(:,1:option%nflowdof)= -ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !22

 
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
 ! print *,'end inter flux'
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_flux.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 0
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#endif
  
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
 !print *,'end jac'
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 ! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
#if 0
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
#endif
#endif

  if (patch%aux%Mphase%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%Mphase%n_zero_rows, &
                          patch%aux%Mphase%zero_rows_local_ghosted,f_up,ierr) 
  endif

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
  print *,'zero rows=', n_zero_rows
  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))
  print *,'zero rows allocated' 
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
print *,'zero rows point 1'
  patch%aux%Mphase%n_zero_rows = n_zero_rows
print *,'zero rows point 2'
  patch%aux%Mphase%zero_rows_local => zero_rows_local
print *,'zero rows point 3'  
  patch%aux%Mphase%zero_rows_local_ghosted => zero_rows_local_ghosted
print *,'zero rows point 4'
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     PETSC_COMM_WORLD,ierr)
  if (flag > 0) patch%aux%Mphase%inactive_cells_exist = .true.

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
 print *,'zero rows', flag
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
  use Patch_module
  use Field_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal :: dcmax, dsmax, max_c, max_S  
  PetscInt :: ierr 

  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  option%dpmax=0.D0
  option%dtmpmax=0.D0 
  option%dcmax=0.D0
  option%dsmax=0.D0
  dcmax=0.D0
  dsmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)

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

  if(option%commsize >1)then
    call MPI_ALLREDUCE(dcmax, max_c,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dsmax, max_s,1, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    dcmax= max_C
    dsmax = max_s
  endif 
  option%dcmax=dcmax
  option%dsmax=dsmax
  !print *, 'Max changes=', option%dpmax,option%dtmpmax, option%dcmax,option%dsmax
end subroutine MphaseMaxChange




subroutine MphaseMaxChangePatch(realization,  max_c, max_s)

  use Realization_module
  use Grid_module
  use Patch_module
  use Field_module
  use Option_module
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: max_s, max_c 


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(option_type), pointer :: option 
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_p(:), yy_p(:),  iphase_loc_p(:),  iphase_old_loc_p(:) 
  PetscInt :: local_id, ghosted_id, n0 
  PetscReal :: cmp
  
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  option => realization%option

  max_c=0.D0
  max_s=0.D0
  
 

  call VecGetArrayF90(field%flow_xx,xx_p, ierr)   
  call VecGetArrayF90(field%flow_yy,yy_p, ierr) 
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p, ierr)
  call VecGetArrayF90(field%iphas_old_loc,iphase_old_loc_p, ierr)

  do local_id =1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif
    n0 = (local_id-1)*option%nflowdof 
    if(int(iphase_loc_p(ghosted_id)) == int(iphase_old_loc_p(ghosted_id)))then
       cmp=dabs(xx_p(n0+3)-yy_p(n0+3))
       if(int(iphase_loc_p(ghosted_id))==1 .or.int(iphase_loc_p(ghosted_id))==2)then
          if(max_c<cmp) max_c = cmp
       endif
       if(int(iphase_loc_p(ghosted_id))==3)then
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
           '"S(l)",' // &
           '"S(g)",' // &
           '"u(l)",'//&
           '"u(g)",'

  do i=1,option%nflowspec
    write(string2,'('',"Xl('',i2,'')"'')') i
    string = trim(string) // trim(string2)
  enddo

  do i=1,option%nflowspec
    write(string2,'('',"Xg('',i2,'')"'')') i
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
  type(mphase_auxvar_type), pointer :: aux_vars(:)
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
 ! print *,'MphaseGetVarFromArray, get pointer'
  if (.not.patch%aux%Mphase%aux_vars_up_to_date) call MphaseUpdateAuxVars(realization)
 ! print *,'MphaseGetVarFromArray, updated'
  aux_vars => patch%aux%Mphase%aux_vars
 ! print *,'MphaseGetVarFromArray, get var'
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
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar)
          case(GAS_ENERGY)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) ! still need implementation
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%den(2)
          case(LIQUID_MOLE_FRACTION)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
          case(LIQUID_ENERGY)
             vec_ptr(local_id) = aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
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
  type(mphase_auxvar_type), pointer :: aux_vars(:)  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  if (.not.patch%aux%Mphase%aux_vars_up_to_date) call MphaseUpdateAuxVars(realization)

  aux_vars => patch%aux%mphase%aux_vars

  select case(ivar)
  case(TEMPERATURE)
    value = aux_vars(ghosted_id)%aux_var_elem(0)%temp
  case(PRESSURE)
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%pres
  case(LIQUID_SATURATION)
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%sat(1)
  case(LIQUID_DENSITY)
      value = aux_vars(ghosted_id)%aux_var_elem(0)%den(1)
  case(GAS_SATURATION)
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%sat(2)
  case(GAS_MOLE_FRACTION)
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(4)
  case(GAS_ENERGY)
      value = aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
  case(GAS_DENSITY) ! still need implementation
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%den(2)
  case(LIQUID_MOLE_FRACTION)
      value  = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2)
  case(LIQUID_ENERGY)
      value = aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
  end select
   
  mphaseGetVarFromArrayAtCell = value
  
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
  !call MphaseAuxDestroy(patch%aux%mphase)

end subroutine MphaseDestroy


#if 0
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

#endif

end module Mphase_module
