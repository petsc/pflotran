module Immis_module
  
  use Immis_Aux_module
  use Global_Aux_module

  implicit none
  
  private 

#include "definitions.h"
  
!#include "include/petscf90.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
!#ifdef USE_PETSC216
!#include "finclude/petscsles.h"
!#endif
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsys.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscerror.h"

! Cutoff parameters
  PetscReal, parameter :: formeps   = 1.D-4
  PetscReal, parameter :: eps = 1.D-5 
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps   = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0  !1D-8
  

  PetscInt, parameter :: jh2o=1, jco2=2

  PetscReal, allocatable, save :: Resold_AR(:,:), Resold_FL(:,:), delx(:,:)
  
  public ImmisResidual,ImmisJacobian, &
         ImmisUpdateFixedAccumulation,ImmisTimeCut,&
         ImmisSetup,ImmisUpdateReason,&
         ImmisMaxChange, ImmisUpdateSolution, &
         ImmisGetTecplotHeader, ImmisInitializeTimestep, &
         ImmisUpdateAuxVars

contains

! ************************************************************************** !
!
! ImmisTimeCut: Resets arrays for time step cut
! author: Chuan Lu
! date: 9/13/08
!
! ************************************************************************** !
subroutine ImmisTimeCut(realization)
 
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

end subroutine ImmisTimeCut

! ************************************************************************** !
!
! ImmisSetup: 
! author: Chuan Lu
! date: 9/13/08
!
! ************************************************************************** !
subroutine ImmisSetup(realization)

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
      call ImmisSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine ImmisSetup

! ************************************************************************** !
!
! ImmisSetupPatch: Creates arrays for auxilliary variables
! author: Chuan Lu
! date: 10/1/08
!
! ************************************************************************** !
subroutine ImmisSetupPatch(realization)

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

  PetscInt :: ghosted_id, iconn, sum_connection, ipara
  type(Immis_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  print *,' ims setup get patch'
  patch%aux%Immis => ImmisAuxCreate()
  
!  option%io_buffer = 'Before Immis can be run, the thc_parameter object ' // &
!                     'must be initialized with the proper variables ' // &
!                     'ImmisAuxCreate() is called anyhwere.'
!  call printErrMsg(option)
  print *,' ims setup get Aux', option%nphase, size(realization%saturation_function_array)     
! immis_parameters create *********************************************
! Sir
  allocate(patch%aux%Immis%Immis_parameter%sir(option%nphase, &
                                  size(realization%saturation_function_array)))
   print *,' ims setup get patch: sir, allocated'                                
  do ipara = 1, size(realization%saturation_function_array)
    patch%aux%Immis%immis_parameter%sir(:,realization%saturation_function_array(ipara)%ptr%id) = &
      realization%saturation_function_array(ipara)%ptr%Sr(:)
  enddo
  print *,' ims setup get patch: sir'
! dencpr  
  allocate(patch%aux%Immis%Immis_parameter%dencpr(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    patch%aux%Immis%Immis_parameter%dencpr(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%rock_density*&
      realization%material_property_array(ipara)%ptr%specific_heat
  enddo
! ckwet
  allocate(patch%aux%Immis%Immis_parameter%ckwet(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    patch%aux%Immis%Immis_parameter%ckwet(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%thermal_conductivity_wet
  enddo
! immis_parameters create_end *****************************************

! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))
  print *,' ims setup get Aux alloc', grid%ngmax
  do ghosted_id = 1, grid%ngmax
    call ImmisAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%Immis%aux_vars => aux_vars
  patch%aux%Immis%num_aux = grid%ngmax
  print *,' ims setup get Aux init'

  allocate(delx(option%nflowdof, grid%ngmax))
  allocate(Resold_AR(grid%nlmax,option%nflowdof))
  allocate(Resold_FL(ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list),option%nflowdof))
  print *,' ims setup allocate app array'
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
  print *,' ims setup get AuxBc alloc', sum_connection
  do iconn = 1, sum_connection
    call ImmisAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%Immis%aux_vars_bc => aux_vars_bc
  patch%aux%Immis%num_aux_bc = sum_connection
  option%numerical_derivatives = PETSC_TRUE

  print *,' ims setup get AuxBc point'
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call ImmisCreateZeroArray(patch,option)

end subroutine ImmisSetupPatch

! ************************************************************************** !
! Immisinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
  function  ImmisInitGuessCheck(realization)
 
  use Realization_module
  use Level_module
  use Patch_module
  use Option_module
  
  PetscInt ::  ImmisInitGuessCheck
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
      ipass= ImmisInitGuessCheckPatch(realization)
      if(ipass<=0)then
        nullify(cur_level)
        exit 
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

   call MPI_Barrier(option%mycomm,ierr)
   if(option%mycommsize >1)then
      call MPI_ALLREDUCE(ipass,ipass0,ONE_INTEGER, MPI_INTEGER,MPI_SUM, &
           option%mycomm,ierr)
      if(ipass0 < option%mycommsize) ipass=-1
   endif
   ImmisInitGuessCheck =ipass
 end function ImmisInitGuessCheck

! ************************************************************************** !
! Immisinitguesscheckpatch: 
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
 subroutine ImmisUpdateReasonPatch(reason,realization)
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
     call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
     call GridVecGetArrayF90(grid,field%flow_yy, yy_p, ierr)

     do n = 1,grid%nlmax
!**** clu-Ignore inactive cells with inactive materials **************
        if (associated(patch%imat)) then
           if (patch%imat(grid%nL2G(n)) <= 0) cycle
        endif
        n0=(n-1)* option%nflowdof
  
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
           if(xx_p(n0 + 3) > 1.D0)then
              re=0; exit
           endif
           if(xx_p(n0 + 3) < 0.)then
              re=0; exit
           endif
     end do
  
    if(re<=0) print *,'Sat out of Region at: ',n,iipha,xx_p(n0+1:n0+3)
    call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call GridVecRestoreArrayF90(grid,field%flow_yy, yy_p, ierr)
    call GridVecRestoreArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr); 

   endif
  
 end subroutine ImmisUpdateReasonPatch


! ************************************************************************** !
!
! ImmisUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Richards problem
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine ImmisUpdateReason(reason, realization)

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
      call ImmisUpdateReasonPatch(re, realization)
        if(re<=0)then
           nullify(cur_level)
           exit 
        endif
        cur_patch => cur_patch%next
     enddo
    cur_level => cur_level%next
 enddo

 call MPI_Barrier(realization%option%mycomm,ierr)
  
  if(realization%option%mycommsize >1)then
     call MPI_ALLREDUCE(re, re0,1, MPI_INTEGER,MPI_SUM, &
          realization%option%mycomm,ierr)
     if(re0<realization%option%mycommsize) re=0
  endif
  reason=re
  
  if(reason<=0 .and. realization%option%myrank ==0) print *,'Sat or Con out of Region', re
end subroutine ImmisUpdateReason

! ************************************************************************** !
! Immisinitguesscheckpatch: 
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
  function  ImmisInitGuessCheckPatch(realization)
   
     use span_wagner_module
     
    use Realization_module
    use Patch_module
    use Field_module
    use Grid_module
    use Option_module
    implicit none
    
    PetscInt :: ImmisInitGuessCheckPatch 
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
    
    call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr)
    
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

    call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
    ImmisInitGuessCheckPatch = ipass
  end function ImmisInitGuessCheckPatch

! ***************************************************************************
!
! ImmisUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Immis problem
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine ImmisUpdateAuxVars(realization)

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
      call ImmisUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine ImmisUpdateAuxVars

! ************************************************************************** !
!
! ImmisUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                        the Immis problem
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
subroutine ImmisUpdateAuxVarsPatch(realization)

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
  type(Immis_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  aux_vars => patch%aux%Immis%aux_vars
  aux_vars_bc => patch%aux%Immis%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    if(.not. associated(realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr))then
       print*, 'error!!! saturation function not allocated', ghosted_id,icap_loc_p(ghosted_id)
    endif
   
    call ImmisAuxVarCompute_NINC(xx_loc_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_var_elem(0), &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       realization%fluid_properties,option)
 ! update global variables
    if( associated(global_aux_vars))then
     
      global_aux_vars(ghosted_id)%pres(:)= aux_vars(ghosted_id)%aux_var_elem(0)%pres -&
               aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
      global_aux_vars(ghosted_id)%temp=aux_vars(ghosted_id)%aux_var_elem(0)%temp
      global_aux_vars(ghosted_id)%sat(:)=aux_vars(ghosted_id)%aux_var_elem(0)%sat(:)
      global_aux_vars(ghosted_id)%den(:)=aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
      global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
                                          * aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
    else
      print *,'Not associated global for IMS'
    endif


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
      case(NEUMANN_BC,ZERO_GRADIENT_BC)
         xxbc(:) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
      end select
      enddo
 
     call ImmisAuxVarCompute_NINC(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0), &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         realization%fluid_properties, option)

     if( associated(global_aux_vars_bc))then
        global_aux_vars_bc(sum_connection)%pres(:)= aux_vars_bc(sum_connection)%aux_var_elem(0)%pres -&
                     aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
        global_aux_vars_bc(sum_connection)%temp=aux_vars_bc(sum_connection)%aux_var_elem(0)%temp
        global_aux_vars_bc(sum_connection)%sat(:)=aux_vars_bc(sum_connection)%aux_var_elem(0)%sat(:)
        !    global_aux_vars(ghosted_id)%sat_store = 
        global_aux_vars_bc(sum_connection)%den(:)=aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
        global_aux_vars_bc(sum_connection)%den_kg = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
                                          * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
  !    global_aux_vars(ghosted_id)%den_kg_store
  !    global_aux_vars(ghosted_id)%mass_balance 
  !    global_aux_vars(ghosted_id)%mass_balance_delta                   
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo


  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  
  patch%aux%Immis%aux_vars_up_to_date = PETSC_TRUE

end subroutine ImmisUpdateAuxVarsPatch

! ************************************************************************** !
!
! ImmisInitializeTimestep: Update data in module prior to time step
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine ImmisInitializeTimestep(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call ImmisUpdateFixedAccumulation(realization)

end subroutine ImmisInitializeTimestep

! ************************************************************************** !
!
! ImmisUpdateSolution: Updates data in module after a successful time step
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine ImmisUpdateSolution(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%flow_xx,realization%field%flow_yy,ierr)   

! make room for hysteric s-Pc-kr

end subroutine ImmisUpdateSolution


! ************************************************************************** !
!
! ImmisUpdateFixedAccumulation: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine ImmisUpdateFixedAccumulation(realization)

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
      call ImmisUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine ImmisUpdateFixedAccumulation

! ************************************************************************** !
!
! ImmisUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine ImmisUpdateFixedAccumPatch(realization)

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
  type(Immis_parameter_type), pointer :: immis_parameter
  type(Immis_auxvar_type), pointer :: aux_vars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  
  call ImmisUpdateAuxVarsPatch(realization) 
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid


  immis_parameter => patch%aux%Immis%immis_parameter
  aux_vars => patch%aux%Immis%aux_vars
    
  call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)

  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call ImmisAccumulation(aux_vars(ghosted_id)%aux_var_elem(0), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,0, accum_p(istart:iend)) 
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)

  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)

#if 0
!  call ImmisNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine ImmisUpdateFixedAccumPatch


! ************************************************************************** !
!
! ImmisAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !  
subroutine ImmisAccumulation(aux_var,por,vol,rock_dencpr,option,iireac,Res)

  use Option_module
  
  implicit none

  type(Immis_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
     
  PetscInt :: ispec, np, iireac
  PetscReal :: porXvol, mol(option%nflowspec), eng
  
 ! if (present(ireac)) iireac=ireac

  porXvol = por*vol
  mol=0.d0; eng=0.D0
  do np = 1, option%nphase
        mol(np) = mol(np) + aux_var%sat(np) * &
             aux_var%den(np)
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
      Res(1:option%nphase)=mol(:)
      Res(option%nflowdof)=eng
  ! endif
  end subroutine ImmisAccumulation

! ************************************************************************** !
!
! ImmisAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !  
subroutine ImmisSourceSink(mmsrc,psrc,tsrc,hsrc,aux_var,isrctype,Res, energy_flag, option)

  use Option_module
  
   use water_eos_module
!   use gas_eos_module  
  use co2eos_module
  use span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use span_wagner_module 
  
  implicit none

  type(Immis_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal mmsrc(option%nflowspec), psrc(option%nphase),tsrc,hsrc 
  PetscInt isrctype
  PetscTruth :: energy_flag
     
  PetscReal :: msrc(option%nflowspec),dw_kg, dw_mol,dddt,dddp
  PetscReal :: enth_src_h2o, enth_src_co2 
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: ukvr, v_darcy, dq, dphi
  PetscInt  :: np, ierr  
  
  Res=0D0
 ! if (present(ireac)) iireac=ireac
      if (energy_flag) then
        Res(option%nflowdof) = Res(option%nflowdof) + hsrc * option%flow_dt   
      endif         
 
   select case(isrctype)
     case(MASS_RATE_SS)
        msrc(:)=mmsrc(:)
             if (msrc(1) > 0.d0) then ! H2O injection
        call wateos_noderiv(tsrc,aux_var%pres,dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        Res(jh2o) = Res( jh2o) + msrc(1) *option%flow_dt
        if (energy_flag) &
             Res(option%nflowdof) = Res(option%nflowdof) + msrc(1)*enth_src_h2o*option%flow_dt
      endif  
    
      if (msrc(2) > 0.d0) then ! CO2 injection
!        call printErrMsg(option,"concentration source not yet implemented in Immis")
      if(option%co2eos == EOS_SPAN_WAGNER)then
         !  span-wagner
          rho = aux_var%den(jco2)*FMWCO2  
          select case(option%itable)  
            case(0,1,2,4,5)
              if( option%itable >=4) then
              call co2_sw_interp(aux_var%pres*1.D-6,&
                  tsrc,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              else
              call co2_span_wagner(aux_var%pres*1.D-6,&
                  tsrc+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              endif 
             case(3) 
              call sw_prop(tsrc,aux_var%pres*1.D-6,rho, &
                     enth_src_co2, eng, fg)
          end select     

         !  units: rho [kg/m^3]; csrc1 [kmol/s]
            enth_src_co2 = enth_src_co2 * FMWCO2
      else if(option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]
            call CO2(tsrc,aux_var%pres, rho,fg, xphi,enth_src_co2)
            enth_src_co2 = enth_src_co2*FMWCO2*option%scale
      else
         call printErrMsg(option,'pflow Immis ERROR: Need specify CO2 EOS')
      endif
              
      Res(jco2) = Res(jco2) + msrc(2)*option%flow_dt
      if (energy_flag) &
         Res(option%nflowdof) = Res(option%nflowdof)+ msrc(2) * enth_src_co2 *option%flow_dt
       endif

     case(-1) ! production well
     ! if node pessure is lower than the given extraction pressure, shut it down
         Dq = psrc(2) ! well parameter, read in input file
                      ! Take the place of 2nd parameter 
        ! Flow term
        do np = 1, option%nphase
          dphi = aux_var%pres - aux_var%pc(np)- psrc(1)
          if (dphi>=0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              v_darcy=0D0
              if (ukvr*Dq>floweps) then
                 v_darcy = Dq * ukvr * dphi
                 Res(np) =Res(np)- v_darcy* aux_var%den(np) 
                 if(energy_flag) Res(option%nflowdof) =Res(option%nflowdof)- v_darcy* aux_var%den(np)*aux_var%h(np)
              endif
           endif
        enddo
       ! print *,'well-prod: ',  aux_var%pres,psrc(1), res
         
    case(1) ! injetion well with constant pressure
         Dq = psrc(2) ! well parameter, read in input file
                      ! Take the place of 2nd parameter 
        ! Flow term
        do np = 1, option%nphase
          dphi = psrc(1) - aux_var%pres - aux_var%pc(np)
          if (dphi>=0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              v_darcy=0D0
              if (ukvr*Dq>floweps) then
                 v_darcy = Dq * ukvr * dphi
                 Res(np) =Res(np)- v_darcy* aux_var%den(np) 
                 if(energy_flag) Res(option%nflowdof) =Res(option%nflowdof)- v_darcy* aux_var%den(np)*aux_var%h(np)
               endif
           endif
        enddo 
    case default
        print *,'Unrecognized Source/Sink condition: ', isrctype 
   end select      
      
 end subroutine ImmisSourceSink


! ************************************************************************** !
!
! ImmisFlux: Computes the internal flux terms for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** ! 
subroutine ImmisFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(Immis_auxvar_elem_type) :: aux_var_up, aux_var_dn
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
!  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
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
        else
           ukvr = aux_var_dn%kvr(np)
           ! if(option%use_isoth == PETSC_FALSE)&
           uh = aux_var_dn%h(np)
        endif
   

        if (ukvr>floweps) then
           v_darcy= Dq * ukvr * dphi
           vv_darcy(np)=v_darcy
           q = v_darcy * area
           fluxm(np)=fluxm(np) + q * density_ave
          ! if(option%use_isoth == PETSC_FALSE)&
            fluxe = fluxe + q*density_ave*uh 
        endif
     endif

#if 0 
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
#endif
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
     Res(1:option%nphase) = fluxm(:) * option%flow_dt
     Res(option%nflowdof) = fluxe * option%flow_dt
 ! end if
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

   end subroutine ImmisFlux

! ************************************************************************** !
!
! ImmisBCFlux: Computes the  boundary flux terms for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine ImmisBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(Immis_auxvar_elem_type) :: aux_var_up, aux_var_dn
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
!  diffdp = por_dn*tor_dn/dd_up*area
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
           v_darcy = aux_vars(MPH_PRESSURE_DOF)
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
        ! uxmol(:)=aux_var_up%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     else
         !if(option%use_isoth == PETSC_FALSE)&
        uh = aux_var_dn%h(np)
         ! uxmol(:)=aux_var_dn%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     endif
    
        fluxm(np) = fluxm(np) + q*density_ave ! *uxmol(ispec)

      !if(option%use_isoth == PETSC_FALSE) &
      fluxe = fluxe + q*density_ave*uh
 !print *,'FLBC', ibndtype(1),np, ukvr, v_darcy, uh, uxmol
   enddo

#if 0 
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
#endif
  ! Conduction term
! if(option%use_isoth == PETSC_FALSE) then
    select case(ibndtype(2))
    case(DIRICHLET_BC, 4)
       Dk =  Dk_dn / dd_up
       cond = Dk*area*(aux_var_up%temp - aux_var_dn%temp) 
       fluxe=fluxe + cond
    end select
! end if

  Res(1:option%nphase)=fluxm(:)* option%flow_dt
  Res(option%nflowdof)=fluxe * option%flow_dt

end subroutine ImmisBCFlux

! ************************************************************************** !
!
! ImmisResidual: Computes the residual equation 
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine ImmisResidual(snes,xx,r,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use grid_module 

  implicit none

  interface
     subroutine samrpetscobjectstateincrease(vec)
       implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       Vec :: vec
     end subroutine samrpetscobjectstateincrease
     
  end interface

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
!  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
 ! check initial guess -----------------------------------------------
  ierr = ImmisInitGuessCheck(realization)
  if(ierr<0)then
    !ierr = PETSC_ERR_ARG_OUTOFRANGE
    if (option%myrank==0) print *,'table out of range: ',ierr
    call SNESSetFunctionDomainError() 
    return
  endif 
  ! end check ---------------------------------------------------------

  ! Communication -----------------------------------------
  ! These 3 must be called before ImmisUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
!  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
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
      call ImmisResidualPatch(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif

end subroutine ImmisResidual

! ************************************************************************** !
!
! ImmisResidualPatch: Computes the residual equation at patch level
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine ImmisResidualPatch(snes,xx,r,realization,ierr)

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
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(Immis_parameter_type), pointer :: immis_parameter
  type(Immis_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscTruth :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  immis_parameter => patch%aux%Immis%immis_parameter
  aux_vars => patch%aux%Immis%aux_vars
  aux_vars_bc => patch%aux%Immis%aux_vars_bc

 ! call ImmisUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%ImmisAux%aux_vars_up_to_date = PETSC_FALSE 

! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%flow_yy,yy_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
!  call GridVecGetArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)
 
 
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
     call ImmisAuxVarCompute_Ninc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0),&
          realization%saturation_function_array(int(icap_loc_p(ng)))%ptr,&
          realization%fluid_properties,option)

     if (option%numerical_derivatives) then
        delx(1,ng) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac * 1.D-3
        delx(2,ng) = xx_loc_p((ng-1)*option%nflowdof+2)*dfac
 
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
         call ImmisAuxVarCompute_Winc(xx_loc_p(istart:iend),delx(:,ng),&
            aux_vars(ng)%aux_var_elem(1:option%nflowdof),&
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
    call ImmisAccumulation(aux_vars(ghosted_id)%aux_var_elem(0),porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
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
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif
      
    psrc(:) = source_sink%flow_condition%pressure%dataset%cur_value(:)
!    qsrc1 = source_sink%flow_condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)
!    hsrc1=0D0
!    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
!    csrc1 = csrc1 / FMWCO2
!    msrc(1)=qsrc1; msrc(2) =csrc1
     msrc(:)= psrc(:)
     msrc(1) =  msrc(1) / FMWH2O
     msrc(2) =  msrc(2) / FMWCO2

     cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      call ImmisSourceSink(msrc,psrc,tsrc1,hsrc1,aux_vars(ghosted_id)%aux_var_elem(0),&
                            source_sink%flow_condition%itype(1),Res,enthalpy_flag, option)
 
      r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o)-Res(jh2o)
      r_p((local_id-1)*option%nflowdof + jco2) = r_p((local_id-1)*option%nflowdof + jco2)-Res(jco2)
      Resold_AR(local_id,jh2o)= Resold_AR(local_id,jh2o) - Res(jh2o)    
      Resold_AR(local_id,jco2)= Resold_AR(local_id,jco2) - Res(jco2)    
      if (enthalpy_flag)then
        r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - Res(option%nflowdof)
        Resold_AR(local_id,option%nflowdof)= Resold_AR(local_id,option%nflowdof) - Res(option%nflowdof)
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
      D_dn = immis_parameter%ckwet(ithrm_dn)

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
!          iphase = int(iphase_loc_p(ghosted_id))
       end select
    enddo

 
    call ImmisAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)

    call ImmisBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(0), &
         aux_vars(ghosted_id)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         immis_parameter%sir(:,icap_dn), &
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
   
      D_up = immis_parameter%ckwet(ithrm_up)
      D_dn = immis_parameter%ckwet(ithrm_dn)

      call ImmisFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
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
  !call GridVecRestoreArrayF90(grid,r, r_p, ierr)


  if (patch%aux%Immis%inactive_cells_exist) then
    do i=1,patch%aux%Immis%n_zero_rows
      r_p(patch%aux%Immis%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_yy, yy_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecRestoreArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
!  call GridVecRestoreArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'Rresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'Rxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
end subroutine ImmisResidualPatch

! ************************************************************************** !
!
! ImmisJacobian: Computes the Jacobian
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine ImmisJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Level_module
  use Grid_module
  use Option_module

  implicit none

  interface
     subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
       
       Mat :: mat
       PetscFortranAddr :: patch
     end subroutine SAMRSetCurrentJacobianPatch
  end interface

  SNES :: snes
  Vec :: xx
  Mat :: A, B, J
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      grid => cur_patch%grid
      ! need to set the current patch in the Jacobian operator
      ! so that entries will be set correctly
      if(associated(grid%structured_grid) .and. &
        (.not.(grid%structured_grid%p_samr_patch.eq.0))) then
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif
      call ImmisJacobianPatch(snes,xx,A,B,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine ImmisJacobian

! ************************************************************************** !
!
! ImmisJacobianPatch: Computes the Jacobian
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine ImmisJacobianPatch(snes,xx,A,B,flag,realization,ierr)

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
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscTruth :: enthalpy_flag
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
  type(Immis_parameter_type), pointer :: immis_parameter
  type(Immis_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  
  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscReal :: psrc(1:realization%option%nphase)
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

  immis_parameter => patch%aux%Immis%immis_parameter
  aux_vars => patch%aux%Immis%aux_vars
  aux_vars_bc => patch%aux%Immis%aux_vars_bc
  
! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
  flag = SAME_NONZERO_PATTERN

#if 0
!  call ImmisNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call GridVecGetArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  call GridVecGetArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
!  call GridVecGetArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)

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
        call ImmisAccumulation(aux_vars(ghosted_id)%aux_var_elem(nvar), &
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
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
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif

    psrc(:) = source_sink%flow_condition%pressure%dataset%cur_value(:)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
 !   hsrc1=0.D0
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)

   ! qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
   ! csrc1 = csrc1 / FMWCO2
      msrc(:)= psrc(:)
      msrc(1) =  msrc(1) / FMWH2O
      msrc(2) =  msrc(2) / FMWCO2
 
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
     do nvar =1, option%nflowdof
       call ImmisSourceSink(msrc,psrc,tsrc1,hsrc1,aux_vars(ghosted_id)%aux_var_elem(nvar),&
                            source_sink%flow_condition%itype(1), Res,enthalpy_flag, option)
      
       ResInc(local_id,jh2o,nvar)=  ResInc(local_id,jh2o,nvar) - Res(jh2o)
       ResInc(local_id,jco2,nvar)=  ResInc(local_id,jco2,nvar) - Res(jco2)
       if (enthalpy_flag) & 
           ResInc(local_id,option%nflowdof,nvar)=&
           ResInc(local_id,option%nflowdof,nvar)- Res(option%nflowdof) 

     enddo 
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
      D_dn = immis_parameter%ckwet(ithrm_dn)

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

 
    call ImmisAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)
    call ImmisAuxVarCompute_Winc(xxbc,delxbc,&
         aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),&
         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties,option)
    
    do nvar=1,option%nflowdof
       call ImmisBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
         aux_vars(ghosted_id)%aux_var_elem(nvar), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         immis_parameter%sir(:,icap_dn), &
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
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
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
      D_up = immis_parameter%ckwet(ithrm_up)
      D_dn = immis_parameter%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
      
      do nvar = 1, option%nflowdof 
         call ImmisFlux(aux_vars(ghosted_id_up)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
            ra(:,nvar)= (Res(:)-ResOld_FL(iconn,:))/delx(nvar,ghosted_id_up)

         call ImmisFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          aux_vars(ghosted_id_dn)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_dn),&
                          tor_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
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
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 0
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#endif
  
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

   
  call GridVecRestoreArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
! call GridVecRestoreArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)
! print *,'end jac'
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

  if (patch%aux%Immis%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%Immis%n_zero_rows, &
                          patch%aux%Immis%zero_rows_local_ghosted,f_up,ierr) 
  endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(option%mycomm,'Rjacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
  endif
end subroutine ImmisJacobianPatch



! ************************************************************************** !
!
! ImmisCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine ImmisCreateZeroArray(patch,option)

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
  patch%aux%Immis%n_zero_rows = n_zero_rows
print *,'zero rows point 2'
  patch%aux%Immis%zero_rows_local => zero_rows_local
print *,'zero rows point 3'  
  patch%aux%Immis%zero_rows_local_ghosted => zero_rows_local_ghosted
print *,'zero rows point 4'
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)
  if (flag > 0) patch%aux%Immis%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
 print *,'zero rows', flag
end subroutine ImmisCreateZeroArray

! ************************************************************************** !
!
! ImmisMaxChange: Computes the maximum change in the solution vector
! author: Chuan Lu
! date: 01/15/08
!
! ************************************************************************** !
subroutine ImmisMaxChange(realization)

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
  PetscReal :: dsmax, max_S  
  PetscInt :: ierr 

  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  option%dpmax=0.D0
  option%dtmpmax=0.D0 
  option%dcmax=0.D0
  option%dsmax=0.D0
  dsmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)
  call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,option%dsmax,ierr)

#if 0
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call ImmisMaxChangePatch(realization, max_s)
      if(dsmax <max_s)  dsmax =max_s
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(option%mycommsize >1)then
    call MPI_ALLREDUCE(dsmax, max_s,1, MPI_DOUBLE_PRECISION,MPI_MAX, option%mycomm,ierr)
    dsmax = max_s
  endif 
  option%dsmax=dsmax
#endif
  !print *, 'Max changes=', option%dpmax,option%dtmpmax, option%dcmax,option%dsmax
end subroutine ImmisMaxChange

! ************************************************************************** !
!
! ImmisGetTecplotHeader: Returns Richards contribution to 
!                               Tecplot file header
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
function ImmisGetTecplotHeader(realization, icolumn)

  use Realization_module
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: ImmisGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-PHASE"'')') icolumn
  else
    write(string2,'('',"PHASE"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-S(l)"'')') icolumn
  else
    write(string2,'('',"S(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-S(g)"'')') icolumn
  else
    write(string2,'('',"S(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-u(l)"'')') icolumn
  else
    write(string2,'('',"u(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-u(g)"'')') icolumn
  else
    write(string2,'('',"u(g)"'')')
  endif
  string = trim(string) // trim(string2)

  ImmisGetTecplotHeader = string

end function ImmisGetTecplotHeader

! ************************************************************************** !
!
! ImmisDestroy: Deallocates variables associated with Immis
! author: Chuan Lu 
! date: 10/14/08
!
! ************************************************************************** !
subroutine ImmisDestroy(patch)

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  !call ImmisAuxDestroy(patch%aux%Immis)

end subroutine ImmisDestroy


#if 0
! ************************************************************************** !
!
! ImmisCheckpointWrite: Writes vecs to checkpoint file
! author: Chuan Lu
! date: 
!
! ************************************************************************** !
subroutine ImmisCheckpointWrite(discretization, viewer)

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call VecView(global_var,viewer,ierr)
  call VecDestroy(global_var,ierr)
  
  
end subroutine ImmisCheckpointWrite


! ************************************************************************** !
!
! ImmisCheckpointRead: Reads vecs from checkpoint file
! author: Chuan Lu 
! date: 
!
! ************************************************************************** !
subroutine ImmisCheckpointRead(discretization,viewer)

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call VecLoadIntoVector(viewer, global_var, ierr)
  call VecDestroy(global_var,ierr)
  
end subroutine ImmisCheckpointRead

#endif

end module Immis_module
