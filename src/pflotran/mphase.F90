module Mphase_module
  
  use Mphase_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
!#include "include/petscf90.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
!#ifdef USE_PETSC216
!#include "finclude/petscsles.h"
!#endif
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscerror.h"

! Cutoff parameters
  PetscReal, parameter :: formeps   = 1.D-4
  PetscReal, parameter :: eps = 1.D-8 
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0  !1D-8
  

  PetscInt, parameter :: jh2o=1, jco2=2

  public MphaseResidual,MphaseJacobian, &
         MphaseUpdateFixedAccumulation,MphaseTimeCut, &
         MphaseSetup,MphaseUpdateReason, &
         MphaseMaxChange,MphaseUpdateSolution, &
         MphaseGetTecplotHeader,MphaseInitializeTimestep, &
         MphaseUpdateAuxVars, init_span_wanger, &
         MphaseSecondaryHeat, MphaseSecondaryHeatJacobian, &
         MphaseComputeMassBalance,MphaseDestroy

contains

! ************************************************************************** !
!
! MphaseTimeCut: Resets arrays for time step cut
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseTimeCut(realization)
 
  use Realization_class
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
! init_span_wanger
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine init_span_wanger(realization)
  use Realization_class
  use co2_span_wagner_module
  use co2_sw_module
  use co2_span_wagner_spline_module

  implicit none
  type(realization_type) :: realization
  PetscMPIInt :: myrank

  if (realization%option%co2eos == EOS_SPAN_WAGNER) then
    select case(realization%option%itable)
       case(0,1,2)
         call initialize_span_wagner(realization%option%itable, &
                                     realization%option%myrank, &
                                     realization%option)
       case(4,5)
         myrank = realization%option%myrank
         call initialize_span_wagner(ZERO_INTEGER,myrank, &
                                     realization%option)
         call initialize_sw_interp(realization%option%itable,myrank)
       case(3)
         call sw_spline_read
       case default
         print *, 'Wrong table option : STOP'
         stop
    end select
  endif
end subroutine init_span_wanger


! ************************************************************************** !
!
! MphaseSetup: 
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseSetup(realization)

  use Realization_class
  use Patch_module
   
  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  call MphaseSetPlotVariables(realization)

end subroutine MphaseSetup

! ************************************************************************** !
!
! MphaseSetupPatch: Creates arrays for auxiliary variables
! author: Chuan Lu
! date: 5/13/08
!
! ************************************************************************** !
subroutine MphaseSetupPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(mphase_type),pointer :: mphase
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection, ipara
  type(Mphase_auxvar_type), pointer :: aux_vars(:)
  type(Mphase_auxvar_type), pointer :: aux_vars_bc(:)
  type(Mphase_auxvar_type), pointer :: aux_vars_ss(:)  
  type(sec_heat_type), pointer :: mphase_sec_heat_vars(:)
  type(coupler_type), pointer :: initial_condition
  PetscReal :: area_per_vol

  option => realization%option
  patch => realization%patch
  grid => patch%grid 

  patch%aux%Mphase => MphaseAuxCreate()
  patch%aux%SC_heat => SecondaryAuxHeatCreate(option)
  mphase => patch%aux%Mphase

  
!  option%io_buffer = 'Before Mphase can be run, the thc_parameter object ' // &
!                     'must be initialized with the proper variables ' // &
!                     'MphaseAuxCreate() is called anywhere.'
!  call printErrMsg(option)  

! mphase_parameters create *********************************************
! Sir
  allocate(mphase%Mphase_parameter%sir(option%nphase, &
                                  size(realization%saturation_function_array)))
  do ipara = 1, size(realization%saturation_function_array)
    mphase%mphase_parameter%sir(:,realization%saturation_function_array(ipara)%ptr%id) = &
      realization%saturation_function_array(ipara)%ptr%Sr(:)
  enddo

! dencpr  
  allocate(mphase%Mphase_parameter%dencpr(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    mphase%mphase_parameter%dencpr(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%rock_density*option%scale*&
      realization%material_property_array(ipara)%ptr%specific_heat
  enddo
! ckwet
  allocate(mphase%Mphase_parameter%ckwet(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    mphase%mphase_parameter%ckwet(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%thermal_conductivity_wet*option%scale
  enddo
  

! mphase_parameters create_end *****************************************

! Create secondary continuum variables - Added by SK 06/28/12

  if (option%use_mc) then
 
    initial_condition => patch%initial_conditions%first
    allocate(mphase_sec_heat_vars(grid%ngmax))
  
    do ghosted_id = 1, grid%ngmax
  
    ! Assuming the same secondary continuum for all regions (need to make it an array)
    ! S. Karra 07/18/12
      call SecondaryContinuumSetProperties( &
        mphase_sec_heat_vars(ghosted_id)%sec_continuum, &
        realization%material_property_array(1)%ptr%secondary_continuum_name, &
        realization%material_property_array(1)%ptr%secondary_continuum_length, &
        realization%material_property_array(1)%ptr%secondary_continuum_matrix_block_size, &
        realization%material_property_array(1)%ptr%secondary_continuum_fracture_spacing, &
        realization%material_property_array(1)%ptr%secondary_continuum_radius, &
        realization%material_property_array(1)%ptr%secondary_continuum_area, &
        option)
        
      mphase_sec_heat_vars(ghosted_id)%ncells = &
        realization%material_property_array(1)%ptr%secondary_continuum_ncells
      mphase_sec_heat_vars(ghosted_id)%aperture = &
        realization%material_property_array(1)%ptr%secondary_continuum_aperture
      mphase_sec_heat_vars(ghosted_id)%epsilon = &
        realization%material_property_array(1)%ptr%secondary_continuum_epsilon
      mphase_sec_heat_vars(ghosted_id)%log_spacing = &
        realization%material_property_array(1)%ptr%secondary_continuum_log_spacing
      mphase_sec_heat_vars(ghosted_id)%outer_spacing = &
        realization%material_property_array(1)%ptr%secondary_continuum_outer_spacing
        

      allocate(mphase_sec_heat_vars(ghosted_id)%area(mphase_sec_heat_vars(ghosted_id)%ncells))
      allocate(mphase_sec_heat_vars(ghosted_id)%vol(mphase_sec_heat_vars(ghosted_id)%ncells))
      allocate(mphase_sec_heat_vars(ghosted_id)%dm_minus(mphase_sec_heat_vars(ghosted_id)%ncells))
      allocate(mphase_sec_heat_vars(ghosted_id)%dm_plus(mphase_sec_heat_vars(ghosted_id)%ncells))
      allocate(mphase_sec_heat_vars(ghosted_id)%sec_continuum% &
               distance(mphase_sec_heat_vars(ghosted_id)%ncells))
    
      call SecondaryContinuumType(&
                              mphase_sec_heat_vars(ghosted_id)%sec_continuum, &
                              mphase_sec_heat_vars(ghosted_id)%ncells, &
                              mphase_sec_heat_vars(ghosted_id)%area, &
                              mphase_sec_heat_vars(ghosted_id)%vol, &
                              mphase_sec_heat_vars(ghosted_id)%dm_minus, &
                              mphase_sec_heat_vars(ghosted_id)%dm_plus, &
                              mphase_sec_heat_vars(ghosted_id)%aperture, &
                              mphase_sec_heat_vars(ghosted_id)%epsilon, &
                              mphase_sec_heat_vars(ghosted_id)%log_spacing, &
                              mphase_sec_heat_vars(ghosted_id)%outer_spacing, &
                              area_per_vol,option)
                                
      mphase_sec_heat_vars(ghosted_id)%interfacial_area = area_per_vol* &
          (1.d0 - mphase_sec_heat_vars(ghosted_id)%epsilon)* &
          realization%material_property_array(1)%ptr% &
          secondary_continuum_area_scaling


    ! Setting the initial values of all secondary node temperatures same as primary node 
    ! temperatures (with initial dirichlet BC only) -- sk 06/26/12
      allocate(mphase_sec_heat_vars(ghosted_id)%sec_temp(mphase_sec_heat_vars(ghosted_id)%ncells))
      
      if (option%set_secondary_init_temp) then
        mphase_sec_heat_vars(ghosted_id)%sec_temp = &
          realization%material_property_array(1)%ptr%secondary_continuum_init_temp
      else
        mphase_sec_heat_vars(ghosted_id)%sec_temp = &
        initial_condition%flow_condition%temperature%dataset%rarray(1)
      endif
          
      mphase_sec_heat_vars(ghosted_id)%sec_temp_update = PETSC_FALSE

    enddo
      
    patch%aux%SC_heat%sec_heat_vars => mphase_sec_heat_vars    
  
  endif

  
  ! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))

  do ghosted_id = 1, grid%ngmax
    call MphaseAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  mphase%aux_vars => aux_vars
  mphase%num_aux = grid%ngmax

  allocate(mphase%delx(option%nflowdof,grid%ngmax))
  allocate(mphase%res_old_AR(grid%nlmax,option%nflowdof))
  allocate(mphase%res_old_FL(ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list),option%nflowdof))

#ifdef YE_FLUX
  allocate(patch%internal_fluxes(3,1,ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list)))
  patch%internal_fluxes = 0.d0
#endif
           
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

  do iconn = 1, sum_connection
    call MphaseAuxVarInit(aux_vars_bc(iconn),option)
  enddo

  mphase%aux_vars_bc => aux_vars_bc
  mphase%num_aux_bc = sum_connection
  
 ! Allocate source /sink  
  source_sink => patch%source_sinks%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    sum_connection = sum_connection + &
                     source_sink%connection_set%num_connections
    source_sink => source_sink%next
  enddo
  allocate(aux_vars_ss(sum_connection))
  do iconn = 1, sum_connection
    call MphaseAuxVarInit(aux_vars_ss(iconn),option)
  enddo
  mphase%aux_vars_ss => aux_vars_ss
  mphase%num_aux_ss = sum_connection
  
  option%numerical_derivatives_flow = PETSC_TRUE

! print *,' mph setup get AuxBc point'
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call MphaseCreateZeroArray(patch,option)

end subroutine MphaseSetupPatch

! ************************************************************************** !
!
! MphaseComputeMassBalance: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine MphaseComputeMassBalance(realization,mass_balance,mass_trapped)

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: mass_trapped(realization%option%nphase)

  type(patch_type), pointer :: cur_patch
  
  mass_balance = 0.d0
  mass_trapped = 0.d0
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseComputeMassBalancePatch(realization,mass_balance,mass_trapped)
    cur_patch => cur_patch%next
  enddo

end subroutine MphaseComputeMassBalance

! ************************************************************************** !
!
! MphaseComputeMassBalancePatch: Initializes mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine MphaseComputeMassBalancePatch(realization,mass_balance,mass_trapped)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
! use Saturation_Function_module
! use Mphase_pckr_module
 
  implicit none
  
  type(realization_type) :: realization
! type(saturation_function_type) :: saturation_function_type

  PetscReal :: mass_balance(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: mass_trapped(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(mphase_auxvar_type), pointer :: mphase_aux_vars(:)
  PetscReal, pointer :: volume_p(:), porosity_loc_p(:), icap_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase
  PetscInt :: ispec_start, ispec_end, ispec
  PetscReal :: pckr_sir(realization%option%nphase)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  mphase_aux_vars => patch%aux%MPhase%aux_vars

  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    ! mass = volume * saturation * density * mole fraction
    do iphase = 1, option%nphase
      do ispec = 1, option%nflowspec
        mass_balance(ispec,iphase) = mass_balance(ispec,iphase) + &
          mphase_aux_vars(ghosted_id)%aux_var_elem(0)%xmol(ispec+(iphase-1)*option%nflowspec)* &
          mphase_aux_vars(ghosted_id)%aux_var_elem(0)%den(iphase)* &
          mphase_aux_vars(ghosted_id)%aux_var_elem(0)%sat(iphase)* &
          porosity_loc_p(ghosted_id)*volume_p(local_id)
      enddo

      pckr_sir(iphase) = &
        realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr%sr(iphase)

      if (iphase == 1 .and. &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%sat(iphase) <= pckr_sir(iphase)) then
        ispec = 1
        mass_trapped(iphase) = mass_trapped(iphase) + &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%xmol(ispec+(iphase-1)*option%nflowspec)* &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%den(iphase)* &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%sat(iphase)* &
        porosity_loc_p(ghosted_id)*volume_p(local_id)
      endif

      if (iphase == 2 .and. &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%sat(iphase) <= pckr_sir(iphase)) then
        ispec = 2
        mass_trapped(iphase) = mass_trapped(iphase) + &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%xmol(ispec+(iphase-1)*option%nflowspec)* &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%den(iphase)* &
        mphase_aux_vars(ghosted_id)%aux_var_elem(0)%sat(iphase)* &
        porosity_loc_p(ghosted_id)*volume_p(local_id)
      endif
    enddo
  enddo

  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p, ierr)
  
end subroutine MphaseComputeMassBalancePatch

! ************************************************************************** !
!
! MphaseZeroMassBalDeltaPatch: Zeros mass balance delta array
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine MphaseZeroMassBalDeltaPatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Mphase%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Mphase%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Mphase%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  
  if (patch%aux%Mphase%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Mphase%num_aux_ss
      global_aux_vars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine MphaseZeroMassBalDeltaPatch

! ************************************************************************** !
!
! MphaseUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine MphaseUpdateMassBalancePatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Mphase%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance = &
      patch%aux%Global%aux_vars(iconn)%mass_balance + &
      patch%aux%Global%aux_vars(iconn)%mass_balance_delta* &
      option%flow_dt
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Mphase%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Mphase%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance = &
        global_aux_vars_bc(iconn)%mass_balance + &
        global_aux_vars_bc(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif
  
  if (patch%aux%Mphase%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Mphase%num_aux_ss
      global_aux_vars_ss(iconn)%mass_balance = &
        global_aux_vars_ss(iconn)%mass_balance + &
        global_aux_vars_ss(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif

end subroutine MphaseUpdateMassBalancePatch

! ************************************************************************** !
! Mphaseinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
function MphaseInitGuessCheck(realization)
 
  use Realization_class
  use Patch_module
  use Option_module
  
  PetscInt ::  MphaseInitGuessCheck
  type(realization_type) :: realization
  type(option_type), pointer:: option
  type(patch_type), pointer :: cur_patch
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr    

  option => realization%option
  ipass = 1
  cur_patch => realization%patch_list%first
  do while(ipass > 0)
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    ipass = MphaseInitGuessCheckPatch(realization)
    cur_patch => cur_patch%next
  enddo

  call MPI_Barrier(option%mycomm,ierr)
  if (option%mycommsize > 1) then
    call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                         option%mycomm,ierr)
    if(ipass0 < option%mycommsize) ipass=-1
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
  use Realization_class
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
  PetscInt :: re0, iipha
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field  
  patch => realization%patch
  grid => patch%grid

  re = 1
 
  if (re > 0) then
    call VecGetArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_yy, yy_p, ierr)
    call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); 
  
    do n = 1,grid%nlmax
!**** clu-Ignore inactive cells with inactive materials **************
      if (associated(patch%imat)) then
        if (patch%imat(grid%nL2G(n)) <= 0) cycle
      endif
      n0 = (n-1)* option%nflowdof
      iipha=int(iphase_loc_p(grid%nL2G(n)))
  
! ******** Too huge change in pressure ****************     
      if (dabs(xx_p(n0 + 1) - yy_p(n0 + 1)) > (1000.0D0 * option%dpmxe)) then
        re = 0; print *,'huge change in p', xx_p(n0 + 1), yy_p(n0 + 1)
        exit
      endif

! ******** Too huge change in temperature ****************
      if (dabs(xx_p(n0 + 2) - yy_p(n0 + 2)) > (10.0D0 * option%dtmpmxe)) then
        re = 0; print *,'huge change in T', xx_p(n0 + 2), yy_p(n0 + 2)
        exit
      endif
 
! ******* Check 0 <= sat/con <= 1 **************************
      select case(iipha)
        case (1)
          if (xx_p(n0 + 3) > 1.0D0) then
            re = 0; exit
          endif
          if (xx_p(n0 + 3) < 0D0) then
            if (xx_p(n0 + 3) > -1D-14) then
              xx_p(n0 + 3) = 0.D0
            else
!             print *,'MPhaseUpdate: ',iipha,n,n0,option%nflowdof,xx_p(n0+3)
              re = 0; exit
            endif          ! clu removed 05/02/2011
          endif
        case (2)
          if (xx_p(n0 + 3) > 1.0D0) then
            re=0; exit
          endif
          if (xx_p(n0 + 3) < 0D-0) then
            if (xx_p(n0 + 3) > -1D-14) then
              xx_p(n0 + 3) = 0.D0
            else  
              re = 0; exit
            endif 
          endif
        case (3)
          if (xx_p(n0 + 3) > 1.D0) then
            re=0; exit
          endif
          if (xx_p(n0 + 3) < 0.) then
            if (xx_p(n0 + 3) > -1D-14) then
              xx_p(n0 + 3) = 0.D0
            else  
              re = 0; exit
            endif  
          endif
      end select
    end do
  
!   if (re <= 0) print *,'Sat or Con out of Region at: ',n,iipha,xx_p(n0+1:n0+3)
    call VecRestoreArrayF90(field%flow_xx, xx_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
    call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr); 

  endif
  ! print *,' update reason', grid%myrank, re,n,grid%nlmax
  reason=re
  
end subroutine MPhaseUpdateReasonPatch


! ************************************************************************** !
!
! MphaseUpdateAuxVars: Updates the auxiliary variables associated with 
!                        the Mphase problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MPhaseUpdateReason(reason, realization)

  use Realization_class
  use Patch_module
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  PetscInt :: reason

  PetscInt :: re, re0
  PetscErrorCode :: ierr

  re = 1
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MPhaseUpdateReasonPatch(re, realization)
    if (re<=0)then
      exit 
    endif
    cur_patch => cur_patch%next
  enddo


  call MPI_Barrier(realization%option%mycomm,ierr)
  
  if (realization%option%mycommsize > 1) then
    call MPI_Allreduce(re,re0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                        realization%option%mycomm,ierr)
    if (re0<realization%option%mycommsize) re=0
  endif
  reason=re
  
  if(reason <= 0 .and. realization%option%myrank == 0) print *,'Sat or Con out of Region', re
end subroutine MPhaseUpdateReason

! ************************************************************************** !
! Mphaseinitguesscheckpatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
  function  MphaseInitGuessCheckPatch(realization)
   
    use co2_span_wagner_module
     
    use Realization_class
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
      
    PetscInt :: local_id, ghosted_id, ipass
    PetscReal, pointer :: xx_p(:)
    PetscErrorCode :: ierr


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
! MphaseUpdateAuxVars: Updates the auxiliary variables associated with 
!                        the Mphase problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseUpdateAuxVars(realization)

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine MphaseUpdateAuxVars

! ************************************************************************** !
!
! MphaseUpdateAuxVarsPatch: Updates the auxiliary variables associated with 
!                        the Mphase problem
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseUpdateAuxVarsPatch(realization)

  use Realization_class
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
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Mphase_auxvar_type), pointer :: aux_vars(:)
  type(Mphase_auxvar_type), pointer :: aux_vars_bc(:) 
  type(Mphase_auxvar_type), pointer :: aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  PetscReal :: xphi, ynacl, mnacl
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  aux_vars => patch%aux%Mphase%aux_vars
  aux_vars_bc => patch%aux%Mphase%aux_vars_bc
  aux_vars_ss => patch%aux%Mphase%aux_vars_ss

  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
  
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
                                  aux_vars(ghosted_id)%aux_var_elem(0),&
                                  global_aux_vars(ghosted_id), &
                                  iphase, &
        realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                                  realization%fluid_properties,option, xphi)
! update global variables
    if (associated(global_aux_vars)) then
      global_aux_vars(ghosted_id)%pres(:) = aux_vars(ghosted_id)%aux_var_elem(0)%pres
      if (iphase == 3) then ! 2-phase
        global_aux_vars(ghosted_id)%pres(1) = aux_vars(ghosted_id)%aux_var_elem(0)%pres - &
               aux_vars(ghosted_id)%aux_var_elem(0)%pc(1)
      endif

!     print *,'UPdate mphase and global vars ', ghosted_id, &
!        aux_vars(ghosted_id)%aux_var_elem(0)%pc(:),aux_vars(ghosted_id)%aux_var_elem(0)%pres, &
!        global_aux_vars(ghosted_id)%pres(:)

      global_aux_vars(ghosted_id)%temp(:) = aux_vars(ghosted_id)%aux_var_elem(0)%temp
      global_aux_vars(ghosted_id)%sat(:) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(:)
      global_aux_vars(ghosted_id)%fugacoeff(1) = xphi
      global_aux_vars(ghosted_id)%den(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
      global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
                                          * aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
      
      mnacl= global_aux_vars(ghosted_id)%m_nacl(1)
      if(global_aux_vars(ghosted_id)%m_nacl(2)>mnacl) mnacl= global_aux_vars(ghosted_id)%m_nacl(2)
      ynacl = mnacl/(1.d3/FMWH2O + mnacl)
      global_aux_vars(ghosted_id)%xmass(1) = (1.d0-ynacl) &
        *aux_vars(ghosted_id)%aux_var_elem(0)%xmol(1) * FMWH2O &
        /((1.d0-ynacl)*aux_vars(ghosted_id)%aux_var_elem(0)%xmol(1) * FMWH2O &
        +aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2) * FMWCO2 &
        +ynacl*aux_vars(ghosted_id)%aux_var_elem(0)%xmol(1)*FMWNACL)
      global_aux_vars(ghosted_id)%xmass(2) = aux_vars(ghosted_id)%aux_var_elem(0)%xmol(3) * FMWH2O &
        /(aux_vars(ghosted_id)%aux_var_elem(0)%xmol(3) * FMWH2O &
        +aux_vars(ghosted_id)%aux_var_elem(0)%xmol(4) * FMWCO2) 
      global_aux_vars(ghosted_id)%reaction_rate_store(:) = global_aux_vars(ghosted_id)%reaction_rate(:)
      global_aux_vars(ghosted_id)%reaction_rate(:) = 0.D0
     !     global_aux_vars(ghosted_id)%mass_balance 
!     global_aux_vars(ghosted_id)%mass_balance_delta                   
    else
      print *,'Not associated global for mph'
    endif
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

! Added the following by Satish Karra 10/05/11
      do idof = 1, option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(HYDROSTATIC_BC,SEEPAGE_BC)
            xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
            if (idof >= MPH_TEMPERATURE_DOF) then
              xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            endif
          case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            iphase = int(iphase_loc_p(ghosted_id))
        end select
      enddo

      select case(boundary_condition%flow_condition%itype(MPH_CONCENTRATION_DOF))
        case(DIRICHLET_BC,SEEPAGE_BC,HYDROSTATIC_BC)
          iphase = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphase = int(iphase_loc_p(ghosted_id))                               
      end select
	  
      call MphaseAuxVarCompute_NINC(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0), &
                          global_aux_vars_bc(sum_connection),iphase, &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         realization%fluid_properties, option, xphi)
    
      if (associated(global_aux_vars_bc)) then
        global_aux_vars_bc(sum_connection)%pres(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%pres -&
                     aux_vars_bc(sum_connection)%aux_var_elem(0)%pc(:)
        global_aux_vars_bc(sum_connection)%temp(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%temp
        global_aux_vars_bc(sum_connection)%sat(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%sat(:)
        !    global_aux_vars(ghosted_id)%sat_store = 
        global_aux_vars_bc(sum_connection)%fugacoeff(1) = xphi
        global_aux_vars_bc(sum_connection)%den(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
        global_aux_vars_bc(sum_connection)%den_kg(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
                                          * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
!       print *,'xxbc ', xxbc, iphasebc, global_aux_vars_bc(sum_connection)%den_kg(:)
        mnacl= global_aux_vars_bc(sum_connection)%m_nacl(1)
        if(global_aux_vars_bc(sum_connection)%m_nacl(2)>mnacl) mnacl = global_aux_vars_bc(sum_connection)%m_nacl(2)
        ynacl =  mnacl/(1.d3/FMWH2O + mnacl)
        global_aux_vars_bc(sum_connection)%xmass(1) = (1.d0-ynacl)&
                              *aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(1) * FMWH2O&
                              /((1.d0-ynacl)*aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(1) * FMWH2O &
                              +aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(2) * FMWCO2 &
                              +ynacl*aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(1)*FMWNACL)
        global_aux_vars_bc(sum_connection)%xmass(2) = aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(3) * FMWH2O&
                              /(aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(3) * FMWH2O&
                              +aux_vars_bc(sum_connection)%aux_var_elem(0)%xmol(4) * FMWCO2) 
 
   
     
  !    global_aux_vars(ghosted_id)%den_kg_store
  !    global_aux_vars(ghosted_id)%mass_balance 
  !    global_aux_vars(ghosted_id)%mass_balance_delta                   
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo


! source/sinks
  source_sink => patch%source_sinks%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      call MphaseAuxVarCopy(aux_vars(ghosted_id)%aux_var_elem(0), &
                              aux_vars_ss(sum_connection)%aux_var_elem(0),option)
      call GlobalAuxVarCopy(global_aux_vars(ghosted_id), &
                            global_aux_vars_ss(sum_connection),option)

    enddo
    source_sink => source_sink%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  patch%aux%Mphase%aux_vars_up_to_date = PETSC_TRUE

end subroutine MphaseUpdateAuxVarsPatch

! ************************************************************************** !
!
! MphaseInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine MphaseInitializeTimestep(realization)

  use Realization_class
  
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

  use Realization_class
  use Field_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch

  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
  
  call VecCopy(realization%field%flow_xx,realization%field%flow_yy,ierr)   
  call VecCopy(realization%field%iphas_loc,realization%field%iphas_old_loc,ierr)
  
  cur_patch => realization%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

! make room for hysteric s-Pc-kr

end subroutine MphaseUpdateSolution

! ************************************************************************** !
!
! MphaseUpdateSolutionPatch: Updates data in module after a successful time 
!                             step 
! author: Satish Karra, LANL
! written based on RichardsUpdateSolutionPatch
! date: 08/23/11
!
! ************************************************************************** !
subroutine MphaseUpdateSolutionPatch(realization)

  use Realization_class
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call MphaseUpdateMassBalancePatch(realization)
  endif

end subroutine MphaseUpdateSolutionPatch

! ************************************************************************** !
!
! MphaseUpdateFixedAccumulation: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !
subroutine MphaseUpdateFixedAccumulation(realization)

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
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

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Secondary_Continuum_Aux_module


  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(mphase_parameter_type), pointer :: mphase_parameter
  type(mphase_auxvar_type), pointer :: aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(sec_heat_type), pointer :: mphase_sec_heat_vars(:)
  
  PetscInt :: ghosted_id, local_id, istart, iend !, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim
    
  call MphaseUpdateAuxVarsPatch(realization)

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  
 
  mphase_parameter => patch%aux%Mphase%mphase_parameter
  aux_vars => patch%aux%Mphase%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  mphase_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

      
  call VecGetArrayF90(field%flow_xx,xx_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%tortuosity_loc,tor_loc_p,ierr)
  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  vol_frac_prim = 1.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    
    if (option%use_mc) then
      vol_frac_prim = mphase_sec_heat_vars(ghosted_id)%epsilon
    endif
    
!    iphase = int(iphase_loc_p(ghosted_id))
!    call MphaseAuxVarCompute_Ninc(xx_p(istart:iend), &
!                       aux_vars(ghosted_id)%aux_var_elem(0), &
!                       iphase, &
!                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
!                       realization%fluid_properties,option)
!    iphase_loc_p(ghosted_id) = iphase
   ! print *, 'MphaseUpdateFixedAccumPatch1'
    !if(.not.associated(aux_vars(ghosted_id))) print *,'no var'
    if(.not.associated(mphase_parameter%dencpr)) print *,'no para'    
    call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(0), &
                              global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              mphase_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,ZERO_INTEGER,vol_frac_prim, &
                              accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%tortuosity_loc,tor_loc_p,ierr)
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
subroutine MphaseAccumulation(aux_var,global_aux_var,por,vol,rock_dencpr, &
                              option,iireac,vol_frac_prim,Res)

  use Option_module
  
  implicit none

  type(mphase_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
  type(global_auxvar_type) :: global_aux_var
        
  PetscInt :: ispec, np, iireac
  PetscReal :: porXvol, mol(option%nflowspec), eng
  PetscReal :: vol_frac_prim
  
 ! if (present(ireac)) iireac=ireac

  porXvol = por*vol
      
  mol=0.d0; eng=0.D0
  do np = 1, option%nphase
    do ispec = 1, option%nflowspec  
      mol(ispec) = mol(ispec) + aux_var%sat(np) * aux_var%den(np) * &
        aux_var%xmol(ispec + (np-1)*option%nflowspec)
    enddo
    eng = eng + aux_var%sat(np) * aux_var%den(np) * aux_var%u(np)
  enddo
  mol = mol * porXvol
 ! if(option%use_isothermal == PETSC_FALSE) &
  eng = eng * porXvol + (1.d0 - por) * vol * rock_dencpr * aux_var%temp 
 
! Reaction terms here
! Note if iireac > 0, then it is the node global index

  if(option%ntrandof > 0)then 
    if (iireac > 0) then
     !H2O
      mol(1) = mol(1) + vol * global_aux_var%reaction_rate_store(1) * &
               option%flow_dt*1D-3 
      !CO2     
      mol(2) = mol(2) + vol * global_aux_var%reaction_rate_store(2) * &
               option%flow_dt*1D-3
    endif
  endif
  
! if (option%use_isothermal)then
!   Res(1:option%nflowdof) = mol(:)
! else
    Res(1:option%nflowdof-1) = vol_frac_prim * mol(:)
    Res(option%nflowdof) = vol_frac_prim * eng
! endif
end subroutine MphaseAccumulation

! ************************************************************************** !
!
! MphaseSourceSink: Computes the source/sink portion for the residual
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !  
subroutine MphaseSourceSink(mmsrc,nsrcpara,psrc,tsrc,hsrc,csrc,aux_var,isrctype,Res, &
                            qsrc_phase,energy_flag,option)

  use Option_module
  use Water_EOS_module
  use EOS_Water_module
!   use Gas_EOS_module  
  use co2eos_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use co2_span_wagner_module
  
  implicit none

  type(mphase_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: psrc(option%nphase),tsrc,hsrc,csrc 
  PetscInt :: isrctype
  PetscInt :: nsrcpara
  PetscBool :: energy_flag
  PetscReal :: qsrc_phase(:) ! volumetric rate of injection/extraction for each phase
     
  PetscReal, allocatable :: msrc(:)
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: enth_src_h2o, enth_src_co2 
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: ukvr, v_darcy, dq, dphi
  PetscReal :: well_status, well_diameter
  PetscReal :: pressure_bh, well_factor, pressure_max, pressure_min
  PetscReal :: well_inj_water, well_inj_co2
  PetscInt  :: np
  PetscInt :: iflag
  PetscErrorCode :: ierr
  
  Res = 0.D0
  allocate(msrc(nsrcpara))
  msrc = mmsrc(1:nsrcpara)

! if (present(ireac)) iireac=ireac
! if (energy_flag) then
!   Res(option%nflowdof) = Res(option%nflowdof) + hsrc * option%flow_dt   
! endif         

  qsrc_phase = 0.d0
  
  select case(isrctype)
    case(MASS_RATE_SS)
      msrc(1) =  msrc(1) / FMWH2O
      msrc(2) =  msrc(2) / FMWCO2
      if (msrc(1) > 0.d0) then ! H2O injection
        call EOSWaterDensityEnthalpy(tsrc,aux_var%pres,dw_kg,dw_mol, &
                                     enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        Res(jh2o) = Res(jh2o) + msrc(1)*(1.d0-csrc)*option%flow_dt
        Res(jco2) = Res(jco2) + msrc(1)*csrc*option%flow_dt
        if (energy_flag) Res(option%nflowdof) = Res(option%nflowdof) + &
          msrc(1)*enth_src_h2o*option%flow_dt
          
!       print *,'soure/sink: ',msrc,csrc,enth_src_h2o,option%flow_dt,option%nflowdof

        ! store volumetric rate for ss_fluid_fluxes()
        qsrc_phase(1) = msrc(1)/dw_mol
      elseif (msrc(1) < 0.d0) then ! H2O extraction
        call EOSWaterDensityEnthalpy(aux_var%temp,aux_var%pres,dw_kg,dw_mol, &
                                     enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        Res(jh2o) = Res(jh2o) + msrc(1)*(1.d0-csrc)*option%flow_dt
        Res(jco2) = Res(jco2) + msrc(1)*csrc*option%flow_dt
        if (energy_flag) Res(option%nflowdof) = Res(option%nflowdof) + &
          msrc(1)*enth_src_h2o*option%flow_dt
          
!       print *,'soure/sink: ',msrc,csrc,enth_src_h2o,option%flow_dt,option%nflowdof

        ! store volumetric rate for ss_fluid_fluxes()
        qsrc_phase(1) = msrc(1)/dw_mol
      endif  
    
      if (msrc(2) > 0.d0) then ! CO2 injection
!       call printErrMsg(option,"concentration source not yet implemented in Mphase")
        if(option%co2eos == EOS_SPAN_WAGNER) then
         !  span-wagner
          rho = aux_var%den(jco2)*FMWCO2  
          select case(option%itable)  
            case(0,1,2,4,5)
              if(option%itable >=4) then
                call co2_sw_interp(aux_var%pres*1.D-6, &
                  tsrc,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              else
                iflag = 1
                call co2_span_wagner(aux_var%pres*1.D-6, &
                  tsrc+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,iflag,option%itable)
              endif 
            case(3) 
              call sw_prop(tsrc,aux_var%pres*1.D-6,rho, &
                     enth_src_co2, eng, fg)
          end select     

         !  units: rho [kg/m^3]; csrc1 [kmol/s]
          enth_src_co2 = enth_src_co2 * FMWCO2
          
          ! store volumetric rate for ss_fluid_fluxes()
          ! qsrc_phase [m^3/sec] = msrc [kmol/sec] / [kg/m^3] * [kg/kmol]  
          qsrc_phase(2) = msrc(2)*rho/FMWCO2

        else if(option%co2eos == EOS_MRK) then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]
          call CO2(tsrc,aux_var%pres, rho,fg, xphi,enth_src_co2)
          qsrc_phase(2) = msrc(2)*rho/FMWCO2
          enth_src_co2 = enth_src_co2*FMWCO2*option%scale
        else
          call printErrMsg(option,'pflow mphase ERROR: Need specify CO2 EOS')
        endif
              
        Res(jco2) = Res(jco2) + msrc(2)*option%flow_dt
        if (energy_flag) Res(option%nflowdof) = Res(option%nflowdof) + msrc(2) * &
          enth_src_co2 *option%flow_dt
        endif

    case(WELL_SS) ! production well
     !if node pessure is lower than the given extraction pressure, shut it down
    ! Flow term
!  well parameter explaination
!   1. well status. 1 injection; -1 production; 0 shut in
!                   2 rate controled injection (same as rate_ss, with max pressure control, not completed yet) 
!                  -2 rate controled production(not implemented for now) 
!
!   2. well factor [m^3],  the effective permeability [m^2/s]
!   3. bottomhole pressure:  [Pa]
!   4. max pressure: [Pa]
!   5. min pressure: [Pa]   
!   6. conc. of water 
!   7. conc. of Co2 (1 - conc. of water)
!   8. well diameter, not used now
!   9. skin factor, not used now

      well_status = msrc(1)
      well_factor = msrc(2)
      pressure_bh = msrc(3)
      pressure_max = msrc(4)
      pressure_min = msrc(5)
      well_inj_water = msrc(6)
      well_inj_co2 = msrc(7)
    
!     if(pressure_min < 0D0) pressure_min = 0D0 !not limited by pressure lower bound   

    ! production well (well status = -1)
      if( dabs(well_status + 1.D0) < 1.D-1) then
        if(aux_var%pres > pressure_min) then
          Dq = well_factor 
          do np = 1, option%nphase
            dphi = aux_var%pres - aux_var%pc(np) - pressure_bh
            if (dphi >= 0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              if(ukvr < 1.e-20) ukvr = 0.D0
              v_darcy = 0.D0
              if (ukvr*Dq > floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = -1.d0*v_darcy
                Res(1) = Res(1) - v_darcy* aux_var%den(np)* &
                  aux_var%xmol((np-1)*option%nflowspec+1)*option%flow_dt
                Res(2) = Res(2) - v_darcy* aux_var%den(np)* &
                  aux_var%xmol((np-1)*option%nflowspec+2)*option%flow_dt
                if(energy_flag) Res(3) = Res(3) - v_darcy * aux_var%den(np)* &
                  aux_var%h(np)*option%flow_dt
              ! print *,'produce: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif 
     !print *,'well-prod: ',  aux_var%pres,psrc(1), res
    ! injection well (well status = 2)
      if( dabs(well_status - 2.D0) < 1.D-1) then

        call EOSWaterDensityEnthalpy(tsrc,aux_var%pres,dw_kg,dw_mol, &
                                     enth_src_h2o,option%scale,ierr)

        Dq = msrc(2) ! well parameter, read in input file
                      ! Take the place of 2nd parameter 
        ! Flow term
        if( aux_var%pres < pressure_max)then  
          do np = 1, option%nphase
            dphi = pressure_bh - aux_var%pres + aux_var%pc(np)
            if (dphi >= 0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              v_darcy = 0.D0
              if (ukvr*Dq>floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = v_darcy
                Res(1) = Res(1) + v_darcy* aux_var%den(np)* &
!                 aux_var%xmol((np-1)*option%nflowspec+1) * option%flow_dt
                  well_inj_water * option%flow_dt
                Res(2) = Res(2) + v_darcy* aux_var%den(np)* &
!                 aux_var%xmol((np-1)*option%nflowspec+2) * option%flow_dt
                  well_inj_co2 * option%flow_dt
!               if(energy_flag) Res(3) = Res(3) + v_darcy*aux_var%den(np)*aux_var%h(np)*option%flow_dt
                if(energy_flag) Res(3) = Res(3) + v_darcy*aux_var%den(np) * &
                  enth_src_h2o*option%flow_dt
                
!               print *,'inject: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif    
    case default
      print *,'Unrecognized Source/Sink condition: ', isrctype 
  end select      
  deallocate(msrc)
      
end subroutine MphaseSourceSink




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
                        option,vv_darcy,vol_frac_prim,Res)
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
  PetscReal :: vv_darcy(:),area,vol_frac_prim
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),fluxe,q, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
#if 0
! This factor 2/3 is multiplied to get bulk perm k=delta^3/12/l, karra 05/14/2013   
  if (option%use_mc) Dq = Dq*2.d0/3.d0*vol_frac_prim
#endif
  diffdp = (por_up*tor_up * por_dn*tor_dn) / &
    (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area*vol_frac_prim 
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy =0.D0 
  
! Flow term
  do np = 1, option%nphase
    
    if (aux_var_up%sat(np) > sir_up(np) .or. aux_var_dn%sat(np) > sir_dn(np)) then
      upweight = dd_dn/(dd_up+dd_dn)
      if (aux_var_up%sat(np) < eps) then 
        upweight = 0.d0
      else if (aux_var_dn%sat(np) < eps) then 
        upweight = 1.d0
      endif
      density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np) 
        
      gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
             (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
             * dist_gravity

      dphi = aux_var_up%pres - aux_var_dn%pres &
             - aux_var_up%pc(np) + aux_var_dn%pc(np) &
             + gravity

      v_darcy = 0.D0
      ukvr = 0.D0
      uh = 0.D0
      uxmol = 0.D0

      ! note uxmol only contains one phase xmol
      ! upstream weighting
      if (dphi >= 0.D0) then
        ukvr = aux_var_up%kvr(np)
        ! if(option%use_isothermal == PETSC_FALSE) &
        uh = aux_var_up%h(np)
        uxmol(1:option%nflowspec) = &
          aux_var_up%xmol((np-1)*option%nflowspec + 1 : np*option%nflowspec)
      else
        ukvr = aux_var_dn%kvr(np)
      ! if(option%use_isothermal == PETSC_FALSE) &
        uh = aux_var_dn%h(np)
        uxmol(1:option%nflowspec) = &
          aux_var_dn%xmol((np-1)*option%nflowspec + 1 : np*option%nflowspec)
      endif

      if (ukvr > floweps) then
        v_darcy = Dq * ukvr * dphi
        vv_darcy(np) = v_darcy
        q = v_darcy * area
        
        do ispec=1, option%nflowspec 
          fluxm(ispec)=fluxm(ispec) + q * density_ave * uxmol(ispec)
        enddo
      ! if(option%use_isothermal == PETSC_FALSE) &
        fluxe = fluxe + q*density_ave*uh 
      endif
    endif

!   Diffusion term   
!   Note : average rule may not be correct

#ifdef PCL

    if ((aux_var_up%sat(np) >= 1.d0) .and. (aux_var_dn%sat(np) >= 1.d0) .or. &
    (aux_var_up%sat(np) <= 0.d0) .and. (aux_var_dn%sat(np) <= 0.d0) &
    ) then
!     single phase
      difff = diffdp * 0.25D0*(aux_var_up%sat(np) + aux_var_dn%sat(np))* &
             (aux_var_up%den(np) + aux_var_dn%den(np))
      do ispec=1, option%nflowspec
        ind = ispec + (np-1)*option%nflowspec
        fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
                (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
!       print *,'mphaseflux1: ',ind,aux_var_up%diff(ind),aux_var_dn%diff(ind)
      enddo
      
    else

!     two-phase
      difff = diffdp * 0.5D0*(aux_var_up%den(np) + aux_var_dn%den(np))
      do ispec=1, option%nflowspec
        ind = ispec + (np-1)*option%nflowspec
        
        if (aux_var_up%xmol(ind) > aux_var_dn%xmol(ind)) then
          upweight = 1.d0
        else
          upweight = 0.d0
        endif
        
        fluxm(ispec) = fluxm(ispec) + difff * &
          0.5D0 * (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
          (upweight*aux_var_up%sat(np)+(1.d0-upweight)*aux_var_dn%sat(np))* &
          (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
!       print *,'mphaseflux1: ',ind,aux_var_up%diff(ind),aux_var_dn%diff(ind)
      enddo
    endif
    
#else

!   print *,'mphaseflux: ',np,aux_var_up%sat(np),aux_var_dn%sat(np),eps, &
!     aux_var_up%den(np),aux_var_dn%den(np),diffdp

!   if ((aux_var_up%sat(np) > eps) .and. (aux_var_dn%sat(np) > eps)) then
      difff = diffdp * 0.25D0*(aux_var_up%sat(np) + aux_var_dn%sat(np))* &
             (aux_var_up%den(np) + aux_var_dn%den(np))
      do ispec=1, option%nflowspec
        ind = ispec + (np-1)*option%nflowspec
        fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (aux_var_up%diff(ind) + aux_var_dn%diff(ind))* &
                (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
!       print *,'mphaseflux1: ',ind,aux_var_up%diff(ind),aux_var_dn%diff(ind)
      enddo
!   endif

#endif

  enddo

! conduction term
  !if(option%use_isothermal == PETSC_FALSE) then     
  Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
  cond = vol_frac_prim * Dk*area * (aux_var_up%temp-aux_var_dn%temp)
  fluxe = fluxe + cond
 ! end if

  !if(option%use_isothermal)then
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
! MphaseBCFlux: Computes boundary flux terms for the residual function
! author: Chuan Lu
! date: 05/12/08
!
! ************************************************************************** !
subroutine MphaseBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,vol_frac_prim,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(mphase_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(:)
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area, vol_frac_prim
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscReal :: Neuman_total_mass_flux, Neuman_mass_flux_spec(option%nflowspec)
  PetscReal :: mol_total_flux(option%nphase)
  PetscInt :: pressure_bc_type

  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  ukvr = 0.d0
  uh = 0.d0

  mol_total_flux = 0.d0
  uxmol = 0.d0

  diffdp = por_dn*tor_dn/dd_up*area*vol_frac_prim

  ! Flow   
  do np = 1, option%nphase
    pressure_bc_type = ibndtype(MPH_PRESSURE_DOF)

    select case(ibndtype(MPH_PRESSURE_DOF))
        ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        Dq = perm_dn / dd_up
#if 0        
! This factor 2/3 is multiplied to get bulk perm k=delta^3/12/l, karra 05/14/2013   
        if (option%use_mc) Dq = Dq*2.d0/3.d0*vol_frac_prim
#endif
        ! Flow term
        ukvr = 0.D0
        v_darcy = 0.D0

!       print *,'Seepage BC: pc/sat ',np, &
!         aux_var_up%pc(np),aux_var_dn%pc(np),aux_var_up%sat(np),aux_var_dn%sat(np)

       if (aux_var_up%sat(np) > sir_dn(np) .or. aux_var_dn%sat(np) > sir_dn(np)) then
          upweight = 1.D0
          if (aux_var_up%sat(np) < eps) then 
            upweight = 0.d0
          else if (aux_var_dn%sat(np) < eps) then 
            upweight = 1.d0
          endif
          density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np)
           
          gravity = (upweight*aux_var_up%den(np) * aux_var_up%avgmw(np) + &
                (1.D0-upweight)*aux_var_dn%den(np) * aux_var_dn%avgmw(np)) &
                * dist_gravity
       
          dphi = aux_var_up%pres - aux_var_dn%pres &
                - aux_var_up%pc(np) + aux_var_dn%pc(np) + gravity

!         print *,'Seepage BC: press ',np,aux_var_up%pres,aux_var_dn%pres,dphi,gravity

          if ((pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC ) .and. np == 2) then
              ! flow in         ! boundary cell is <= pref
            if (dphi > 0.d0) then
              dphi = 0.d0
            endif
          endif

          if (dphi >= 0.D0) then
            ukvr = aux_var_up%kvr(np)
          else
            ukvr = aux_var_dn%kvr(np)
          endif
     
          if (ukvr*Dq > floweps) then
            v_darcy = Dq * ukvr * dphi
          endif

!         print *,'Seepage: vD/kvr/dphi ',np,V_darcy,ukvr,Dq,dphi
        endif

        q = v_darcy * area 
        vv_darcy(np) = v_darcy

        uh = 0.D0
        uxmol = 0.D0
        mol_total_flux(np) = q*density_ave  
        if (v_darcy >= 0.D0) then
!     if(option%use_isothermal == PETSC_FALSE) &
          uh = aux_var_up%h(np)
          uxmol(:) = aux_var_up%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
        else
!     if(option%use_isothermal == PETSC_FALSE) &
          uh = aux_var_dn%h(np)
          uxmol(:) = aux_var_dn%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
        endif
 
     
      case(NEUMANN_BC)
        v_darcy = 0.D0
        if (dabs(aux_vars(1)) > floweps) then
          Neuman_total_mass_flux = aux_vars(MPH_PRESSURE_DOF)
          if (v_darcy > 0.d0) then 
            density_ave = aux_var_up%den(np)
          else 
            density_ave = aux_var_dn%den(np)
          endif
          if (np == 1) then
            Neuman_mass_flux_spec(np) = &
            Neuman_total_mass_flux * (1.D0-aux_vars(MPH_CONCENTRATION_DOF))
            uxmol(1) = 1.D0; uxmol(2)=0.D0
            mol_total_flux(np) = Neuman_mass_flux_spec(np)/FMWH2O
            uh = aux_var_dn%h(np)
          else
            Neuman_mass_flux_spec(np) = &
            Neuman_total_mass_flux * aux_vars(MPH_CONCENTRATION_DOF)
            uxmol(1) = 0.D0; uxmol(2) = 1.D0
            mol_total_flux(np) = Neuman_mass_flux_spec(np)/FMWCO2
            uh = aux_var_dn%h(np)
          endif
          vv_darcy(np) = mol_total_flux(np)/density_ave
        endif

      case(ZERO_GRADIENT_BC)

    end select
     
     
    
    do ispec=1, option%nflowspec 
      fluxm(ispec) = fluxm(ispec) + mol_total_flux(np)*uxmol(ispec)
    enddo
      !if(option%use_isothermal == PETSC_FALSE) &
    fluxe = fluxe + mol_total_flux(np)*uh
!   print *,'FLBC', ibndtype(1),np, ukvr, v_darcy, uh, uxmol, mol_total_flux
  enddo
  
! Diffusion term   
  select case(ibndtype(MPH_CONCENTRATION_DOF))
    case(DIRICHLET_BC) 
  ! if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
     !diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
      do np = 1, option%nphase
        if (aux_var_up%sat(np)>eps .and. aux_var_dn%sat(np) > eps) then
          diff = diffdp * 0.25D0*(aux_var_up%sat(np)+aux_var_dn%sat(np))* &
                    (aux_var_up%den(np)+aux_var_up%den(np))
          do ispec = 1, option%nflowspec
            fluxm(ispec) = fluxm(ispec) + &
              diff * aux_var_dn%diff((np-1)* option%nflowspec+ispec)* &
              (aux_var_up%xmol((np-1)* option%nflowspec+ispec) &
              -aux_var_dn%xmol((np-1)* option%nflowspec+ispec))
          enddo
        endif         
      enddo
     
  end select

! Conduction term
! if(option%use_isothermal == PETSC_FALSE) then
    select case(ibndtype(MPH_TEMPERATURE_DOF))
      case(DIRICHLET_BC)
        Dk =  Dk_dn / dd_up
        cond = vol_frac_prim * Dk*area*(aux_var_up%temp - aux_var_dn%temp)
        fluxe = fluxe + cond
      case(NEUMANN_BC)
        fluxe = fluxe + aux_vars(MPH_TEMPERATURE_DOF)*area*option%scale
	  ! aux_vars(MPH_TEMPERATURE_DOF) stores heat flux, 1.d-6 is to convert
	  ! from W to MW, Added by Satish Karra, LANL 10/05/11
      case(ZERO_GRADIENT_BC)
      ! No change in fluxe
    end select
!   print *, fluxe, aux_vars
! end if

  Res(1:option%nflowspec) = fluxm(:) * option%flow_dt
  Res(option%nflowdof) = fluxe * option%flow_dt

end subroutine MphaseBCFlux

! ************************************************************************** !
!
! MphaseResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine MphaseResidual(snes,xx,r,realization,ierr)

  use Realization_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module 

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
  type(patch_type), pointer :: cur_patch
  PetscInt :: ichange, i  

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
    call SNESSetFunctionDomainError(snes,ierr) 
    return
  endif 
  ! end check ---------------------------------------------------------


  ! Variable switching-------------------------------------------------
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseVarSwitchPatch(xx, realization, ZERO_INTEGER, ichange)
    call MPI_Allreduce(ichange,i,ONE_INTEGER_MPI,MPIU_INTEGER, &
                        MPI_MIN,option%mycomm,ierr)
    ichange = i 
    if (ichange < 0) then
      call SNESSetFunctionDomainError(snes,ierr) 
      return
    endif
    cur_patch => cur_patch%next
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
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseResidualPatch(snes,xx,r,realization,ierr)
    cur_patch => cur_patch%next
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

  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  
  use EOS_Water_module
  use Gas_EOS_module  
  use co2eos_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use co2_span_wagner_module

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
  PetscReal :: k1, k2, z1, z2, xg, vmco2, vmh2o, sg, sgg
  PetscReal :: xmol(realization%option%nphase*realization%option%nflowspec),&
               satu(realization%option%nphase)
  PetscReal :: yh2o_in_co2 = 1.d-2
! PetscReal :: yh2o_in_co2 = 0.d0
  PetscReal :: wat_sat_x, co2_sat_x
  PetscReal :: lngamco2, m_na, m_cl, m_nacl, Qkco2, mco2, xco2eq, temp
! PetscReal :: xla,co2_poyn
  PetscInt :: local_id, ghosted_id, dof_offset
  PetscInt :: iflag
  PetscInt :: idum
  PetscReal :: min_value
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars(:)

  patch => realization%patch  
  grid => patch%grid
  option => realization%option
  field => realization%field
  global_aux_vars => patch%aux%Global%aux_vars

#if 0
  option%force_newton_iteration = PETSC_FALSE
  ! checking for negative saturation/mole fraction
  call VecStrideMin(xx,TWO_INTEGER,idum,min_value,ierr)
  if (min_value < 0.d0) then
    write(option%io_buffer,*) 'Warning: saturation or mole fraction negative at cell ', &
      idum, min_value 
    call printMsg(option)
    option%force_newton_iteration = PETSC_TRUE
  endif
#endif
    
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
    iipha=int(iphase_loc_p(ghosted_id))
    p = xx_p(dof_offset+1)
    t = xx_p(dof_offset+2)
    den(1:option%nphase) = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(1:option%nphase)
    select case(iipha) 
      case(1) ! liquid
        xmol(2) = xx_p(dof_offset+3)
        xmol(1) = 1.D0 - xmol(2)
        satu(1) = 1.D0; satu(2) = 0.D0
      case(2) ! gas
        xmol(4) = xx_p(dof_offset+3)
        xmol(3) = 1.D0 - xmol(4)
        satu(1) = 0.d0; satu(2) = 1.D0
      case(3) ! two-phase
        satu(2) = xx_p(dof_offset+3) 
        satu(1) = 1.D0 - satu(2)
        xmol(3) = yh2o_in_co2; xmol(4) = 1.D0-xmol(3)
    end select

! Pure CO2 phase properties ------------------------------------------    
    p2 = p
!   p2 = p*xmol(4)
    if(p2 >= 5.d4)then
      if(option%co2eos == EOS_SPAN_WAGNER)then
        select case(option%itable)  
          case(0,1,2,4,5)
            if(option%itable >=4) then
              call co2_sw_interp(p2*1.D-6,t,dg,dddt,dddp,fg, &
                  dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            else
              iflag = 1
              call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg, &
                  dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag,option%itable)
              if (iflag < 1) then
                ichange = -1
                return
              endif
            endif
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
! Span-Wagner EOS with Bi-Cubic Spline interpolation
          case(3) 
            call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
        end select     
      elseif(option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
        call CO2( t,p2, dg,fg, xphi, hg)
        dg = dg / FMWCO2
        hg = hg * FMWCO2 *option%scale
      endif
    else      
      call ideal_gaseos_noderiv(p2,t,option%scale,dg,hg,ug)
      fg = p2
    endif
   
    xphi = fg/p2
    call EOSWaterSaturationPressure(t, sat_pressure, ierr)
    sat_pressure = sat_pressure/1.D5
  
    m_na=option%m_nacl; m_cl=m_na; m_nacl = m_na 
    if(associated(realization%reaction))then
      if (associated(realization%reaction%species_idx)) then
        if (realization%reaction%species_idx%na_ion_id /= 0 .and. &
          realization%reaction%species_idx%cl_ion_id /= 0) then
          m_na = global_aux_vars(ghosted_id)%m_nacl(1)
          m_cl = global_aux_vars(ghosted_id)%m_nacl(2)
          m_nacl = m_na
          if (m_cl > m_na) m_nacl = m_cl
        endif
      endif  
    endif

    call Henry_duan_sun(t,p2*1.D-5,henry,xphi,lngamco2, &
      m_na,m_cl,sat_pressure)

    Qkco2 = henry*xphi ! QkCO2 = xphi * exp(-mu0) / gamma

    sat_pressure = sat_pressure * 1.D5
    mco2 = (p - sat_pressure)*1.D-5 * Qkco2 ! molality CO2, y * P = P - Psat(T)

    xco2eq = mco2/(1.D3/fmwh2o + mco2 + m_nacl) ! mole fraction CO2

    henry = 1.D8 / FMWH2O / henry / xphi !note: henry = H/phi
    wat_sat_x = sat_pressure/p ! X_w^sc
    co2_sat_x = (1.D0-wat_sat_x)/(henry/p-wat_sat_x)*henry/p  ! xmol(4) = xmol(2)*henry/p

!     tmp = 1.D0-tmp ! approximate form

    select case(icri) ! icri = 0 in call statement
      case(0)
      select case(iipha)     
        case(1) ! liquid
          xmol(4) = xmol(2)*henry/p   

!         print *,'phase chg: ',xmol(2),xco2eq,mco2,m_nacl,p,t

!         if(xmol(4)+ wat_sat_x > 1.05d0) then
!         if(xmol(2) > xco2eq) then
          if(xmol(2) >= xco2eq) then

        !   Rachford-Rice initial guess: 1=H2O, 2=CO2
            k1 = wat_sat_x !sat_pressure*1.D5/p
            k2 = henry/p
            z1 = xmol(1); z2 = xmol(2)
            xg = ((1.d0-k2)*z2+(1.d0-k1)*z1)/((1.d0-k2)*(1.d0-k1)*(z1+z2))
            vmco2 = 1.d0/dg
            vmh2o = 1.D0 /den(1)   ! FMWH2O/0.9d3
            
       !    calculate initial guess for sg
            sg = vmco2*xg/(vmco2*xg+vmh2o*(1.d0-xg))
            if(sg>1.D-4) then          
              write(*,'('' Liq -> 2ph '',''rank='',i6,'' n='',i8,'' p='',1pe10.4, &
       &      '' T='',1pe10.4,'' Xl='',1pe11.4, &
       &      '' Xco2eq='',1pe11.4,'' sg='',1pe11.4)') &
              option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3),xco2eq,sg

              iphase_loc_p(ghosted_id) = 3 ! Liq -> 2ph
        
!             write(*,'(''Rachford-Rice: '','' z1, z2='', &
!       &     1p2e12.4,'' xeq='',1pe12.4,'' xg='',1pe12.4, &
!       &     '' sg='',1pe12.4)') z1,z2,xco2eq,xg,sg
        
!             write(*,'(''Rachford-Rice: '',''K1,2='',1p2e12.4,'' z1,2='', &
!       &     1p2e12.4,'' xeq='',1pe12.4,'' xg='',1pe12.4, &
!       &     '' sg='',1pe12.4,'' dg='',1p2e12.4)') &
!             k1,k2,z1,z2,xco2eq,xg,sg,den(1),dg

!             sgg = den(1)*(z2-xco2eq)/(den(1)*(z2-xco2eq) - &
!               dg*(z2-(1.d0-wat_sat_x)))
!             write(*,'(''Rachford-Rice: sg = '',1p2e12.4)') sgg,sg
              
              xx_p(dof_offset+3) = sg   
              ichange = 1
            endif
          endif

        case(2) ! gas
          
!         if (xmol(3) > wat_sat_x * 1.05d0) then
          if (xmol(3) > wat_sat_x * 1.001d0) then

!           print *,'gas -> 2ph: ',xmol(3),wat_sat_x,xco2eq,sat_pressure
          
!         if (xmol(3) > (1.d0+1.d-6)*tmp .and. iipha==2)then
            write(*,'('' Gas -> 2ph '',''rank='',i6,'' n='',i8, &
       &  '' p= '',1pe10.4,'' T= '',1pe10.4,'' Xg= '',1pe11.4,'' Ps/P='', 1pe11.4)') &
            option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3),wat_sat_x

            iphase_loc_p(ghosted_id) = 3 ! Gas -> 2ph
           xx_p(dof_offset+3) = 1.D0-formeps
!            xx_p(dof_offset+3) = 1.D0
            ichange = 1
          endif

        case(3) ! two-phase
          tmp = wat_sat_x
          
!         xmol(2)= (1.D0-tmp)/(Henry/p-tmp) ! solve: x1+x2=1, y1+y2=1, y1=k1*x1, y2=k2*x2
!         xmol(2)= p*(1.D0-tmp)/Henry ! approximate form
!         xmol(1)= 1.D0-xmol(2)
!         xmol(3)= xmol(1)*tmp
!         xmol(4)= 1.D0-xmol(3)

          xmol(1) = 1.D0 - xco2eq
          xmol(2) = xco2eq
          xmol(3) = tmp
          xmol(4) = 1.D0-tmp          

!         if(satu(2) >= 1.D0) then
!         if(satu(2) >= 1.01D0) then
          if(satu(2) >= 1.001D0) then

            write(*,'('' 2ph -> Gas '',''rank='',i6,'' n='',i8, &
       &  '' p='',1pe10.4,'' T='',1pe10.4,'' sg='',1pe11.4)') &
            option%myrank,local_id,xx_p(dof_offset+1:dof_offset+3)

            iphase_loc_p(ghosted_id) = 2 ! 2ph -> Gas
!           xx_p(dof_offset+3) = 1.D0 - 1.D-8
            xx_p(dof_offset+3) = 1.D0
            ichange = 1
            
          else if(satu(2) <= 0.D0) then
          
            write(*,'('' 2ph -> Liq '',''rank= '',i6,'' n='',i8,'' p='',1pe10.4, &
      &     '' T='',1pe10.4,'' sg ='',1pe11.4,'' Xco2eq='',1pe11.4)')  &
            option%myrank,local_id, xx_p(dof_offset+1:dof_offset+3),xmol(2)

            iphase_loc_p(ghosted_id) = 1 ! 2ph -> Liq
            ichange = 1
            tmp = xmol(2) * 0.99
            xx_p(dof_offset+3) = tmp
          endif

      end select

    case(1)
    
      select case(iipha)     
        case(2)   
          if (xmol(3) > wat_sat_x) then
!           write(*,'(''** Gas -> 2ph '',i8,1p10e12.4)') local_id, &
!           xx_p(dof_offset+1:dof_offset+3)
          endif
        case(1) 
          xmol(4) = xmol(2)*henry/p 
          if (xmol(4) >co2_sat_x ) then
!           write(*,'(''** Liq -> 2ph '',i8,1p10e12.4)') local_id, &
!           xx_p(dof_offset+1:dof_offset+3),xmol(4), co2_sat_x
          endif
        case(3) 
          if(satu(2) > 1.D0 .and. iipha == 3) then
!            write(*,'(''** 2ph -> Gas '',i8,1p10e12.4)') local_id, &
!            xx_p(dof_offset+1:dof_offset+3)
          endif
          if(satu(2) <= 0.D0 .and. iipha == 3) then
!           write(*,'(''** 2ph -> Liq '',i8,1p10e12.4)') local_id, &
!           xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
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
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  
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
  PetscReal :: psrc(1:realization%option%nphase)
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(mphase_type), pointer :: mphase
  type(mphase_parameter_type), pointer :: mphase_parameter
  type(mphase_auxvar_type), pointer :: aux_vars(:)
  type(mphase_auxvar_type), pointer :: aux_vars_bc(:)
  type(mphase_auxvar_type), pointer :: aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal, pointer :: msrc(:)

  type(sec_heat_type), pointer :: mphase_sec_heat_vars(:)

  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: nsrcpara
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: vol_frac_prim
  
  ! secondary continuum variables
  PetscReal :: sec_dencpr
  PetscReal :: area_prim_sec
  PetscReal :: res_sec_heat
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  mphase => patch%aux%Mphase
  mphase_parameter => mphase%mphase_parameter
  aux_vars => mphase%aux_vars
  aux_vars_bc => mphase%aux_vars_bc
  aux_vars_ss => mphase%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

  mphase_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

 ! call MphaseUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MphaseAux%aux_vars_up_to_date = PETSC_FALSE 
 
  if (option%compute_mass_balance_new) then
    call MphaseZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
 
! call VecGetArrayF90(field%flow_yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
 

  vol_frac_prim = 1.d0
  
#if 1
  ! Pertubations for aux terms --------------------------------
  do ng = 1, grid%ngmax
    if (grid%nG2L(ng) < 0) cycle
    if (associated(patch%imat)) then
      if (patch%imat(ng) <= 0) cycle
    endif
        
    istart = (ng-1) * option%nflowdof + 1; iend = istart - 1 + option%nflowdof
    iphase = int(iphase_loc_p(ng))
    ghosted_id = ng
    call MphaseAuxVarCompute_Ninc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0), &
      global_aux_vars(ng), iphase, &
      realization%saturation_function_array(int(icap_loc_p(ng)))%ptr, &
      realization%fluid_properties,option,xphi)

#if 1
    if (associated(global_aux_vars)) then
       global_aux_vars(ghosted_id)%pres(:) = aux_vars(ghosted_id)%aux_var_elem(0)%pres - &
               aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
       global_aux_vars(ghosted_id)%temp(:) = aux_vars(ghosted_id)%aux_var_elem(0)%temp
       global_aux_vars(ghosted_id)%sat(:) = aux_vars(ghosted_id)%aux_var_elem(0)%sat(:)
       global_aux_vars(ghosted_id)%fugacoeff(1) = xphi
       global_aux_vars(ghosted_id)%den(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
       global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
                                          * aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
!       global_aux_vars(ghosted_id)%reaction_rate(:)=0D0
!      print *,'UPdate mphase and gloable vars', ghosted_id, global_aux_vars(ghosted_id) %m_nacl(:), & 
!      global_aux_vars(ghosted_id)%pres(:)
!      global_aux_vars(ghosted_id)%mass_balance 
!      global_aux_vars(ghosted_id)%mass_balance_delta                   
    else
       print *,'Not associated global for mph'
    endif
#endif

    if (option%numerical_derivatives_flow) then
      mphase%delx(1,ng) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac !* 1.D-3
      mphase%delx(2,ng) = xx_loc_p((ng-1)*option%nflowdof+2)*dfac

!     print *,'mphase_delx: ',dfac,mphase%delx(1,ng),mphase%delx(2,ng),mphase%delx(3,ng), &
!     xx_loc_p((ng-1)*option%nflowdof+1),xx_loc_p((ng-1)*option%nflowdof+2),xx_loc_p((ng-1)*option%nflowdof+3)

      select case (iphase)
        case (1)
          if (xx_loc_p((ng-1)*option%nflowdof+3) < 5D-5) then
            mphase%delx(3,ng) =  dfac*xx_loc_p((ng-1)*option%nflowdof+3)
          else
            mphase%delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
          endif
          if (mphase%delx(3,ng) < 1D-8 .and. mphase%delx(3,ng) >= 0.D0) mphase%delx(3,ng) = 1D-8
          if (mphase%delx(3,ng) > -1D-8 .and. mphase%delx(3,ng) < 0.D0) mphase%delx(3,ng) = -1D-8
        case(2)  
          if (xx_loc_p((ng-1)*option%nflowdof+3) < 0.9995) then
            mphase%delx(3,ng) =  dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
          else
            mphase%delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
          endif
          if (mphase%delx(3,ng) < 1D-8 .and. mphase%delx(3,ng) >= 0.D0) mphase%delx(3,ng) = 1D-8
          if (mphase%delx(3,ng) > -1D-8 .and. mphase%delx(3,ng) < 0.D0) mphase%delx(3,ng) = -1D-8
        case(3)
          if (xx_loc_p((ng-1)*option%nflowdof+3) <= 0.9) then
            mphase%delx(3,ng) = dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
          else
            mphase%delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
          endif
           
          if (mphase%delx(3,ng) < 1D-12 .and. mphase%delx(3,ng) >= 0.D0) mphase%delx(3,ng) = 1D-12
          if (mphase%delx(3,ng) > -1D-12 .and. mphase%delx(3,ng) < 0.D0) mphase%delx(3,ng) = -1D-12
        
          if ((mphase%delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3)) > 1.D0) then
            mphase%delx(3,ng) = (1.D0-xx_loc_p((ng-1)*option%nflowdof+3))*1D-6
          endif
          if ((mphase%delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3)) < 0.D0) then
            mphase%delx(3,ng) = xx_loc_p((ng-1)*option%nflowdof+3)*1D-6
          endif
      end select
      call MphaseAuxVarCompute_Winc(xx_loc_p(istart:iend),mphase%delx(:,ng),&
            aux_vars(ng)%aux_var_elem(1:option%nflowdof),global_aux_vars(ng),iphase,&
            realization%saturation_function_array(int(icap_loc_p(ng)))%ptr,&
            realization%fluid_properties,option)
    endif
  enddo
! print *,'mphase resi patch: end numerical increments'
#endif

  mphase%res_old_AR=0.D0; mphase%res_old_FL=0.D0; r_p = 0.d0

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
    
    if (option%use_mc) then
      vol_frac_prim = mphase_sec_heat_vars(ghosted_id)%epsilon
    endif
    
    call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(0), &
                            global_aux_vars(ghosted_id), &
                            porosity_loc_p(ghosted_id), &
                            volume_p(local_id), &
                            mphase_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                            option,ONE_INTEGER,vol_frac_prim,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
 !   print *,'REs, acm: ', res
    mphase%res_old_AR(local_id, :)= Res(1:option%nflowdof)
  enddo
#endif


! ================== Secondary continuum heat source terms =====================
#if 1
  if (option%use_mc) then
  ! Secondary continuum contribution (Added by SK 06/26/2012)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
    
      sec_dencpr = mphase_parameter%dencpr(int(ithrm_loc_p(ghosted_id))) ! secondary rho*c_p same as primary for now

      if (option%sec_vars_update) then
        call MphaseSecHeatAuxVarCompute(mphase_sec_heat_vars(ghosted_id), &
                        aux_vars(ghosted_id)%aux_var_elem(0), &
                        global_aux_vars(ghosted_id), &
                        mphase_parameter%ckwet(int(ithrm_loc_p(ghosted_id))), &
                        sec_dencpr, &
                        option)
      endif       
    
      call MphaseSecondaryHeat(mphase_sec_heat_vars(ghosted_id), &
                        aux_vars(ghosted_id)%aux_var_elem(0), &
                        global_aux_vars(ghosted_id), &
                        mphase_parameter%ckwet(int(ithrm_loc_p(ghosted_id))), &
                        sec_dencpr, &
                        option,res_sec_heat) 
      r_p(iend) = r_p(iend) - res_sec_heat*option%flow_dt*volume_p(local_id)

    enddo   
    option%sec_vars_update = PETSC_FALSE
  endif
#endif

! ============== end secondary continuum heat source ===========================

#if 1
  ! Source/sink terms -------------------------------------
! print *, 'Mphase residual patch 2' 
  source_sink => patch%source_sinks%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    !print *, 'RES s/s begin'
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
    enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif

    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif
!   qsrc1 = source_sink%flow_condition%pressure%dataset%rarray(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)
    
!   print *,'src/sink: ',tsrc1,csrc1,hsrc1,psrc
    
!   hsrc1=0D0
!   qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
!   csrc1 = csrc1 / FMWCO2
!   msrc(1)=qsrc1; msrc(2) =csrc1
!geh begin change
!geh remove
!geh    msrc(:)= psrc(:)
!clu add
  select case(source_sink%flow_condition%itype(1))
    case(MASS_RATE_SS)
      msrc => source_sink%flow_condition%rate%dataset%rarray
      nsrcpara= 2
    case(WELL_SS)
      msrc => source_sink%flow_condition%well%dataset%rarray
      nsrcpara = 7 + option%nflowspec 
     
!    print *,'src/sink: ',nsrcpara,msrc
    case default
      print *, 'mphase mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
      stop  
  end select

!clu end change

    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call MphaseSourceSink(msrc,nsrcpara, psrc,tsrc1,hsrc1,csrc1, &
                            aux_vars(ghosted_id)%aux_var_elem(0),&
                            source_sink%flow_condition%itype(1),Res, &
                    ! fluid flux [m^3/sec] = Res [kmol/mol] / den [kmol/m^3]
                            patch%ss_fluid_fluxes(:,sum_connection), &
                            enthalpy_flag,option)

  ! included by SK, 08/23/11 to print mass fluxes at source/sink						
      if (option%compute_mass_balance_new) then
        global_aux_vars_ss(sum_connection)%mass_balance_delta(:,1) = &
          global_aux_vars_ss(sum_connection)%mass_balance_delta(:,1) - &
          Res(:)/option%flow_dt
      endif
  
      r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o)-Res(jh2o)
      r_p((local_id-1)*option%nflowdof + jco2) = r_p((local_id-1)*option%nflowdof + jco2)-Res(jco2)
      mphase%res_old_AR(local_id,jh2o) = mphase%res_old_AR(local_id,jh2o) - Res(jh2o)    
      mphase%res_old_AR(local_id,jco2) = mphase%res_old_AR(local_id,jco2) - Res(jco2)    
      if (enthalpy_flag) then
        r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - &
          Res(option%nflowdof)
        mphase%res_old_AR(local_id,option%nflowdof) = &
          mphase%res_old_AR(local_id,option%nflowdof) - Res(option%nflowdof)
      endif 
  !   else if (qsrc1 < 0.d0) then ! withdrawal
  !   endif
    enddo
    source_sink => source_sink%next
  enddo
#endif

#if 1
 ! print *, 'Mphase residual patch 3' 
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
      D_dn = mphase_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  
! Then need fill up increments for BCs
      do idof =1, option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(HYDROSTATIC_BC,SEEPAGE_BC)
            xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
            if(idof>=MPH_TEMPERATURE_DOF)then
              xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            endif 
          case(ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            iphase = int(iphase_loc_p(ghosted_id))
          case(NEUMANN_BC)  
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            iphase = int(iphase_loc_p(ghosted_id))
        end select
      enddo

      select case(boundary_condition%flow_condition%itype(MPH_CONCENTRATION_DOF))
        case(DIRICHLET_BC,SEEPAGE_BC,HYDROSTATIC_BC)
          iphase = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphase=int(iphase_loc_p(ghosted_id))                               
      end select
 
      call MphaseAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
            global_aux_vars_bc(sum_connection), iphase,&
            realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
            realization%fluid_properties, option, xphi)

#if 1
      if( associated(global_aux_vars_bc))then
        global_aux_vars_bc(sum_connection)%pres(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%pres -&
                     aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
        global_aux_vars_bc(sum_connection)%temp(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%temp
        global_aux_vars_bc(sum_connection)%sat(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%sat(:)
      !    global_aux_vars(ghosted_id)%sat_store = 
        global_aux_vars_bc(sum_connection)%fugacoeff(1) = xphi
        global_aux_vars_bc(sum_connection)%den(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
        global_aux_vars_bc(sum_connection)%den_kg = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
                                          * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
  !     global_aux_vars(ghosted_id)%den_kg_store
  !     global_aux_vars(ghosted_id)%mass_balance 
  !     global_aux_vars(ghosted_id)%mass_balance_delta                   
      endif
#endif


      call MphaseBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(0), &
         aux_vars(ghosted_id)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tor_loc_p(ghosted_id), &
         mphase_parameter%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,vol_frac_prim,Res)

      patch%boundary_velocities(:,sum_connection) = v_darcy(:)

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)

      mphase%res_old_AR(local_id,1:option%nflowdof) = &
           mphase%res_old_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
   !  print *, 'REs BC: ',r_p(istart:iend)

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_aux_vars_bc(sum_connection)%mass_balance_delta(:,1) = &
          global_aux_vars_bc(sum_connection)%mass_balance_delta(:,1) &
            - Res(:)/option%flow_dt 
      endif

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
                         dot_product(option%gravity, &
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
   
      D_up = mphase_parameter%ckwet(ithrm_up)
      D_dn = mphase_parameter%ckwet(ithrm_dn)


      call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
          tor_loc_p(ghosted_id_up),mphase_parameter%sir(:,icap_up), &
          dd_up,perm_up,D_up, &
          aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
          tor_loc_p(ghosted_id_dn),mphase_parameter%sir(:,icap_dn), &
          dd_dn,perm_dn,D_dn, &
          cur_connection_set%area(iconn),distance_gravity, &
          upweight,option,v_darcy,vol_frac_prim,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy(:)
      mphase%res_old_FL(sum_connection,1:option%nflowdof)= Res(1:option%nflowdof)
      
      if (local_id_up > 0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn > 0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

#ifdef YE_FLUX
      patch%internal_fluxes(1:option%nflowdof,1,sum_connection) = &
                                                     Res(1:option%nflowdof)
#endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif

! adjust residual to R/dt
  select case (option%idt_switch) 
  case(1) 
    r_p(:) = r_p(:)/option%flow_dt
  case(-1)
    if(option%flow_dt > 1.D0) r_p(:) = r_p(:)/option%flow_dt
  end select
  
  do local_id = 1, grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    endif

    istart = 1 + (local_id-1)*option%nflowdof
!   if(volume_p(local_id) > 1.D0) &    ! karra added 05/06/2013
    r_p (istart:istart+option%nflowdof-1) = &
        r_p(istart:istart+option%nflowdof-1)/volume_p(local_id)
    if(r_p(istart) > 1E20 .or. r_p(istart) < -1E20) print *,'mphase residual: ', r_p (istart:istart+2)
  enddo

! print *,'finished rp vol scale'
  if(option%use_isothermal) then
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


  if (mphase%inactive_cells_exist) then
    do i=1,mphase%n_zero_rows
      r_p(mphase%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr)
! call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'Rresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'MPHxx.out',viewer,ierr)
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

  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  
  flag = SAME_NONZERO_PATTERN
  call MatGetType(A,mat_type,ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr)
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseJacobianPatch(snes,xx,J,J,flag,realization,ierr)
    cur_patch => cur_patch%next
  enddo

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'MPHjacobian.out', &
                              viewer,ierr)
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)    
    call MatNorm(J,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)    
    call MatNorm(J,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
  endif

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
  use Realization_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i, j
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
! PetscReal :: max_dev
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: natural_id_up,natural_id_dn
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
               Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof,&
           realization%option%nflowdof)
  PetscReal :: dummy_real_array(2)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(mphase_type), pointer :: mphase
  type(mphase_parameter_type), pointer :: mphase_parameter
  type(mphase_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  
  type(sec_heat_type), pointer :: sec_heat_vars(:)
  
  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal, pointer :: msrc(:)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: iphasebc                
  PetscInt :: nsrcpara  
   
  PetscViewer :: viewer
  Vec :: debug_vec
  PetscReal :: vol_frac_prim

  ! secondary continuum variables
  PetscReal :: area_prim_sec
  PetscReal :: jac_sec_heat
  
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

  mphase => patch%aux%mphase
  mphase_parameter => mphase%mphase_parameter
  aux_vars => mphase%aux_vars
  aux_vars_bc => mphase%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  
  sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  
! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
#if 0
!  call MphaseNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  ResInc = 0.D0
  vol_frac_prim = 1.d0
 
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
     
    if (option%use_mc) then
      vol_frac_prim = sec_heat_vars(ghosted_id)%epsilon
    endif
     
    do nvar =1, option%nflowdof
      call MphaseAccumulation(aux_vars(ghosted_id)%aux_var_elem(nvar), &
             global_aux_vars(ghosted_id), &
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             mphase_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
             option,ONE_INTEGER,vol_frac_prim,res) 
      ResInc(local_id,:,nvar) = ResInc(local_id,:,nvar) + Res(:)
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

    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
 !   hsrc1=0.D0
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)

   ! qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
   ! csrc1 = csrc1 / FMWCO2
!geh begin change
!geh remove
!geh      msrc(:)= psrc(:)
!clu add
     select case(source_sink%flow_condition%itype(1))
     case(MASS_RATE_SS)
       msrc => source_sink%flow_condition%rate%dataset%rarray
       nsrcpara= 2
     case(WELL_SS)
       msrc => source_sink%flow_condition%well%dataset%rarray
       nsrcpara = 7 + option%nflowspec 
     case default
       print *, 'mphase mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
       stop  
     end select

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
        
       do nvar = 1, option%nflowdof
         call MphaseSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,csrc1, &
                               aux_vars(ghosted_id)%aux_var_elem(nvar), &
                               source_sink%flow_condition%itype(1),Res, &
                               dummy_real_array,enthalpy_flag,option)
      
         ResInc(local_id,jh2o,nvar) =  ResInc(local_id,jh2o,nvar) - Res(jh2o)
         ResInc(local_id,jco2,nvar) =  ResInc(local_id,jco2,nvar) - Res(jco2)
         if (enthalpy_flag) & 
           ResInc(local_id,option%nflowdof,nvar) = &
           ResInc(local_id,option%nflowdof,nvar) - Res(option%nflowdof) 

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
      D_dn = mphase_parameter%ckwet(ithrm_dn)
      
 
      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))

! Then need fill up increments for BCs
      delxbc = 0.D0;
      do idof = 1, option%nflowdof   
        select case(boundary_condition%flow_condition%itype(idof))
        case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof) = 0.D0
        case(HYDROSTATIC_BC,SEEPAGE_BC)
          xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
          if (idof >= MPH_TEMPERATURE_DOF) then
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            delxbc(idof) = mphase%delx(idof,ghosted_id)
          endif 
        case(ZERO_GRADIENT_BC,NEUMANN_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          !iphasebc = int(iphase_loc_p(ghosted_id))
          delxbc(idof) = mphase%delx(idof,ghosted_id)
        end select
      enddo
    !print *,'BC:',boundary_condition%flow_condition%itype, xxbc, delxbc


      select case(boundary_condition%flow_condition%itype(MPH_CONCENTRATION_DOF))
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
      case(ZERO_GRADIENT_BC)
        iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select
      if (boundary_condition%flow_condition%itype(MPH_PRESSURE_DOF) /= NEUMANN_BC) then
        call MphaseAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0), &
           global_aux_vars_bc(sum_connection),iphasebc,&
           realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
           realization%fluid_properties,option)
        call MphaseAuxVarCompute_Winc(xxbc,delxbc,&
           aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),&
           global_aux_vars_bc(sum_connection),iphasebc, &
           realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
           realization%fluid_properties,option)
    
        do nvar=1,option%nflowdof
          call MphaseBCFlux(boundary_condition%flow_condition%itype, &
            boundary_condition%flow_aux_real_var(:,iconn), &
            aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
            aux_vars(ghosted_id)%aux_var_elem(nvar), &
            porosity_loc_p(ghosted_id), &
            tor_loc_p(ghosted_id), &
            mphase_parameter%sir(:,icap_dn), &
            cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
            cur_connection_set%area(iconn), &
            distance_gravity,option, &
            vv_darcy,vol_frac_prim,Res)
          ResInc(local_id,1:option%nflowdof,nvar) = &
              ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
        enddo
      endif
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

    ra = 0.D0
!   max_dev = 0.D0
    do neq = 1, option%nflowdof
      do nvar = 1, option%nflowdof
        ra(neq,nvar) = (ResInc(local_id,neq,nvar) - mphase%res_old_AR(local_id,neq)) &
                     /mphase%delx(nvar,ghosted_id)

!       print *,'finite dif: ',neq,nvar,ra(neq,nvar),ResInc(local_id,neq,nvar),mphase%res_old_AR(local_id,neq), &
!       mphase%delx(nvar,ghosted_id)

!       if (max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
      enddo
    enddo
   
    select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) = ra(1:option%nflowdof,1:option%nflowdof)/option%flow_dt
      case(-1)
        if(option%flow_dt > 1) ra(1:option%nflowdof,1:option%nflowdof) = ra(1:option%nflowdof,1:option%nflowdof)/option%flow_dt
    end select

    Jup = ra(1:option%nflowdof,1:option%nflowdof)

    if (option%use_mc) then

      call MphaseSecondaryHeatJacobian(sec_heat_vars(ghosted_id), &
                        mphase_parameter%ckwet(int(ithrm_loc_p(ghosted_id))), &
                        mphase_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                        option,jac_sec_heat)
 ! sk - option%flow_dt cancels out with option%flow_dt in the denominator for the term below                                      
      Jup(option%nflowdof,2) = Jup(option%nflowdof,2) - &
                               jac_sec_heat*volume_p(local_id) 
    endif

!   if (volume_p(local_id) > 1.D0) &    ! karra added 05/06/2013
      Jup = Jup / volume_p(local_id)

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
                         dot_product(option%gravity, &
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
    
      iphas_up = int(iphase_loc_p(ghosted_id_up))
      iphas_dn = int(iphase_loc_p(ghosted_id_dn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      D_up = mphase_parameter%ckwet(ithrm_up)
      D_dn = mphase_parameter%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))

      
      do nvar = 1, option%nflowdof 
        call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(nvar), &
                         porosity_loc_p(ghosted_id_up), &
                         tor_loc_p(ghosted_id_up), &
                         mphase_parameter%sir(:,icap_up), &
                         dd_up,perm_up,D_up, &
                         aux_vars(ghosted_id_dn)%aux_var_elem(0), &
                         porosity_loc_p(ghosted_id_dn), &
                         tor_loc_p(ghosted_id_dn), &
                         mphase_parameter%sir(:,icap_dn), &
                         dd_dn,perm_dn,D_dn, &
                         cur_connection_set%area(iconn), &
                         distance_gravity, &
                         upweight, option, vv_darcy, vol_frac_prim, Res)
                         
        ra(:,nvar) = (Res(:)-mphase%res_old_FL(iconn,:))/ &
                     mphase%delx(nvar,ghosted_id_up)

        call MphaseFlux(aux_vars(ghosted_id_up)%aux_var_elem(0), &
                         porosity_loc_p(ghosted_id_up), &
                         tor_loc_p(ghosted_id_up), &
                         mphase_parameter%sir(:,icap_up), &
                         dd_up,perm_up,D_up, &
                         aux_vars(ghosted_id_dn)%aux_var_elem(nvar), &
                         porosity_loc_p(ghosted_id_dn),&
                         tor_loc_p(ghosted_id_dn), &
                         mphase_parameter%sir(:,icap_dn), &
                         dd_dn,perm_dn,D_dn, &
                         cur_connection_set%area(iconn),distance_gravity, &
                         upweight, option, vv_darcy, vol_frac_prim, Res)
      
        ra(:,nvar+option%nflowdof) = (Res(:)-mphase%res_old_FL(iconn,:)) &
                                      /mphase%delx(nvar,ghosted_id_dn)
      enddo

      select case(option%idt_switch)
      case(1)
        ra = ra / option%flow_dt
      case(-1)  
        if(option%flow_dt>1)  ra = ra / option%flow_dt
      end select
    
      if (local_id_up > 0) then
        voltemp=1.D0
!       if (volume_p(local_id_up) > 1.D0)then   ! karra added 05/06/2013
          voltemp = 1.D0/volume_p(local_id_up)
!       endif
        Jup(:,1:option%nflowdof)= ra(:,1:option%nflowdof)*voltemp !11
        Jdn(:,1:option%nflowdof)= &
          ra(:,1 + option%nflowdof:2 * option%nflowdof)*voltemp !12

        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
        voltemp = 1.D0
!       if (volume_p(local_id_dn) > 1.D0) then   ! karra added 05/06/2013
          voltemp = 1.D0/volume_p(local_id_dn)
!       endif
        Jup(:,1:option%nflowdof) = -ra(:,1:option%nflowdof)*voltemp !21
        Jdn(:,1:option%nflowdof) = -ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !22

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
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tor_loc_p, ierr)
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
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
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

  if (mphase%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,mphase%n_zero_rows, &
                          mphase%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif

end subroutine MphaseJacobianPatch



! ************************************************************************** !
!
! MphaseCreateZeroArray: Computes the zeroed rows for inactive grid cells
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
! print *,'zero rows=', n_zero_rows
  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))
! print *,'zero rows allocated' 
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
!print *,'zero rows point 1'
  patch%aux%Mphase%n_zero_rows = n_zero_rows
!print *,'zero rows point 2'
  patch%aux%Mphase%zero_rows_local => zero_rows_local
!print *,'zero rows point 3'  
  patch%aux%Mphase%zero_rows_local_ghosted => zero_rows_local_ghosted
!print *,'zero rows point 4'
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)
  if (flag > 0) patch%aux%Mphase%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
! print *,'zero rows', flag
end subroutine MphaseCreateZeroArray

! ************************************************************************** !
!
! MphaseMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/15/08
!
! ************************************************************************** !
subroutine MphaseMaxChange(realization)

  use Realization_class
  use Patch_module
  use Field_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  PetscReal :: dcmax, dsmax, max_c, max_S  
  PetscErrorCode :: ierr 

  option => realization%option
  field => realization%field

  option%dpmax=0.D0
  option%dtmpmax=0.D0 
  option%dcmax=0.D0
  option%dsmax=0.D0
  dcmax=0.D0
  dsmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call MphaseMaxChangePatch(realization, max_c, max_s)
    if(dcmax <max_c)  dcmax =max_c
    if(dsmax <max_s)  dsmax =max_s
    cur_patch => cur_patch%next
  enddo

  if(option%mycommsize >1)then
    call MPI_Allreduce(dcmax,max_c,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    call MPI_Allreduce(dsmax,max_s,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    dcmax= max_C
    dsmax = max_s
  endif 
  option%dcmax=dcmax
  option%dsmax=dsmax
  !print *, 'Max changes=', option%dpmax,option%dtmpmax, option%dcmax,option%dsmax
end subroutine MphaseMaxChange




subroutine MphaseMaxChangePatch(realization,  max_c, max_s)

  use Realization_class
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
! MphaseGetTecplotHeader: Returns Mphase contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
function MphaseGetTecplotHeader(realization,icolumn)

  use Realization_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: MphaseGetTecplotHeader
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
    write(string2,'('',"'',i2,''-P(l) [Pa]"'')') icolumn
  else
    write(string2,'('',"P(l) [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P(g) [Pa]"'')') icolumn
  else
    write(string2,'('',"P(g) [Pa]"'')')
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
    write(string2,'('',"'',i2,''-d(l)"'')') icolumn
  else
    write(string2,'('',"d(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-d(g)"'')') icolumn
  else
    write(string2,'('',"d(g)"'')')
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
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-vis(l)"'')') icolumn
  else
    write(string2,'('',"vis(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-vis(g)"'')') icolumn
  else
    write(string2,'('',"vis(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-kvr(l)"'')') icolumn
  else
    write(string2,'('',"kvr(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-kvr(g)"'')') icolumn
  else
    write(string2,'('',"kvr(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i1,'')"'')') icolumn, i
    else
      write(string2,'('',"Xl('',i1,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xg('',i1,'')"'')') icolumn, i
    else
      write(string2,'('',"Xg('',i1,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-PHASE"'')') icolumn
  else
    write(string2,'('',"PHASE"'')')
  endif
  string = trim(string) // trim(string2)
  
  MphaseGetTecplotHeader = string

end function MphaseGetTecplotHeader

! ************************************************************************** !
!
! MphaseSetPlotVariables: Adds variables to be printed to list
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine MphaseSetPlotVariables(realization)
  
  use Realization_class
  use Output_Aux_module
  use Variables_module

  implicit none

  type(realization_type) :: realization
  type(output_variable_type) :: output_variable
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => realization%output_option%output_variable_list
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)
  
  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)
  
  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)

  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)

  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)

  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)

  name = 'Liquid Energy'
  units = 'kJ/mol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)

  name = 'Gas Energy'
  units = 'kJ/mol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)

  name = 'Liquid Viscosity'
  units = 'Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_VISCOSITY)

  name = 'Gas Viscosity'
  units = 'Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_VISCOSITY)

  name = 'Liquid Mobility'
  units = '1/Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOBILITY)

  name = 'Gas Mobility'
  units = '1/Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOBILITY)

  name = 'Liquid Mole Fraction H2O'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION,ONE_INTEGER)

  name = 'Liquid Mole Fraction CO2'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION,TWO_INTEGER)

  name = 'Gas Mole Fraction H2O'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION,ONE_INTEGER)

  name = 'Gas Mole Fraction CO2'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION,TWO_INTEGER)

  name = 'Phase'
  units = ''
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               PHASE)

end subroutine MphaseSetPlotVariables

! ************************************************************************** !
!
! MphaseSecondaryHeat: Calculates the source term contribution due to secondary
! continuum in the primary continuum residual 
! author: Satish Karra, LANL
! date: 06/26/12
!
! ************************************************************************** !
subroutine MphaseSecondaryHeat(sec_heat_vars,aux_var,global_aux_var, &
                            therm_conductivity,dencpr, &
                            option,res_heat)
                            
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module

  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(mphase_auxvar_elem_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  PetscReal :: temp_current_N
  PetscReal :: res_heat

  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
! temp_primary_node = global_aux_var%temp(1)
  temp_primary_node = aux_var%temp

  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! We only need the temperature at the outer-most node (closest to primary node)
  temp_current_N = rhs(ngcells)/coeff_diag(ngcells)
  
  ! Calculate the coupling term
  res_heat = area_fm*therm_conductivity*(temp_current_N - temp_primary_node)/ &
             dm_plus(ngcells)
                          
end subroutine MphaseSecondaryHeat


! ************************************************************************** !
!
! MphaseSecondaryHeatJacobian: Calculates the source term jacobian contribution 
! due to secondary continuum in the primary continuum residual 
! author: Satish Karra, LANL
! date: 06/6/12
!
! ************************************************************************** !
subroutine MphaseSecondaryHeatJacobian(sec_heat_vars, &
                                    therm_conductivity, &
                                    dencpr, &
                                    option,jac_heat)
                                    
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module

  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: m
  PetscReal :: Dtemp_N_Dtemp_prim
  PetscReal :: jac_heat
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  area_fm = sec_heat_vars%interfacial_area
  dm_minus = sec_heat_vars%dm_minus
 
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                                        
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    ! We do not have to calculate rhs terms
  enddo

  ! We need the temperature derivative at the outer-most node (closest to primary node)
  Dtemp_N_Dtemp_prim = 1.d0/coeff_diag(ngcells)*alpha*area(ngcells)/ &
                       (dm_plus(ngcells)*vol(ngcells))
  
  ! Calculate the jacobian term
  jac_heat = area_fm*therm_conductivity*(Dtemp_N_Dtemp_prim - 1.d0)/ &
             dm_plus(ngcells)
                            
              
end subroutine MphaseSecondaryHeatJacobian

! ************************************************************************** !
!
! MphaseDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MphaseDestroy(realization)

  use Realization_class

  implicit none
  
  type(realization_type) :: realization
  
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
  
  call VecLoad(global_var, viewer, ierr)
  call VecDestroy(global_var,ierr)
  ! solid volume fraction
  if (mphase_option%rk > 0.d0) then
    call VecLoad(mphase_field%phis, viewer, ierr)
  endif  
  
end subroutine MphaseCheckpointRead

#endif

end module Mphase_module
