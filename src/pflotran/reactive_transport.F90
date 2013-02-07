module Reactive_Transport_module

  use Transport_module
  use Reaction_module

  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: RTTimeCut, &
            RTSetup, &
            RTMaxChange, &
            RTUpdateSolution, &
            RTResidual, &
            RTJacobian, &
            RTInitializeTimestep, &
            RTUpdateAuxVars, &
            RTComputeMassBalance, &
            RTDestroy, &
            RTUpdateTransportCoefs, &
            RTUpdateRHSCoefs, &
            RTCalculateRHS_t0, &
            RTCalculateRHS_t1, &
            RTCalculateTransportMatrix, &
            RTReact, &
            RTCheckUpdate, &
            RTJumpStartKineticSorption, &
            RTCheckpointKineticSorption, &
            RTExplicitAdvection
  
contains

! ************************************************************************** !
!
! RTTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTTimeCut(realization)
 
  use Realization_module
  use Field_module
  use Global_module
 
  implicit none
  
  type(realization_type) :: realization
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr

  field => realization%field
 
  ! copy previous solution back to current solution
  call VecCopy(field%tran_yy,field%tran_xx,ierr)
  
  ! set densities and saturations to t
  if (realization%option%nflowdof > 0) then
    call GlobalUpdateDenAndSat(realization,realization%option%tran_weight_t0)
  endif
  
  call RTInitializeTimestep(realization)  
  ! note: RTUpdateTransportCoefs() is called within RTInitializeTimestep()
  ! geh - not any longer, tran coefs should always be evaluated at time k+1  
  
  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalUpdateDenAndSat(realization,realization%option%tran_weight_t1)
  endif

  call RTUpdateTransportCoefs(realization)
 
end subroutine RTTimeCut

! ************************************************************************** !
!
! RTSetup: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RTSetup(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTSetupPatch(realization)
  call RTSetPlotVariables(realization)

end subroutine RTSetup

! ************************************************************************** !
!
! RTSetupPatch: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RTSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Constraint_module
  use Fluid_module
  use Material_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(fluid_property_type), pointer :: cur_fluid_property
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(coupler_type), pointer :: initial_condition
  type(tran_constraint_type), pointer :: sec_tran_constraint

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: iphase
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction
  sec_tran_constraint => realization%sec_transport_constraint

  patch%aux%RT => RTAuxCreate(option)
  patch%aux%SC_RT => SecondaryAuxRTCreate(option)
  patch%aux%RT%rt_parameter%ncomp = reaction%ncomp
  patch%aux%RT%rt_parameter%naqcomp = reaction%naqcomp
  patch%aux%RT%rt_parameter%offset_aqueous = reaction%offset_aqueous
  patch%aux%RT%rt_parameter%nimcomp = reaction%nimcomp
  patch%aux%RT%rt_parameter%offset_immobile = reaction%offset_immobile
  if (reaction%ncollcomp > 0) then
    patch%aux%RT%rt_parameter%ncoll = reaction%ncoll
    patch%aux%RT%rt_parameter%offset_colloid  = reaction%offset_colloid 
    patch%aux%RT%rt_parameter%ncollcomp = reaction%ncollcomp
    patch%aux%RT%rt_parameter%offset_collcomp = reaction%offset_collcomp
    allocate(patch%aux%RT%rt_parameter%pri_spec_to_coll_spec(reaction%naqcomp))
    patch%aux%RT%rt_parameter%pri_spec_to_coll_spec = &
      reaction%pri_spec_to_coll_spec
    allocate(patch%aux%RT%rt_parameter%coll_spec_to_pri_spec(reaction%ncollcomp))
    patch%aux%RT%rt_parameter%coll_spec_to_pri_spec = &
      reaction%coll_spec_to_pri_spec
  endif
  if (reaction%nimcomp > 0) then
    patch%aux%RT%rt_parameter%nimcomp = reaction%nimcomp
    patch%aux%RT%rt_parameter%offset_immobile = reaction%offset_immobile
  endif
 
!============== Create secondary continuum variables - SK 2/5/13 ===============

  
  if (option%use_mc) then
    initial_condition => patch%initial_conditions%first
    allocate(rt_sec_transport_vars(grid%ngmax))  
    do ghosted_id = 1, grid%ngmax
    ! Assuming the same secondary continuum type for all regions
      call SecondaryRTAuxVarInit(realization%material_property_array(1)%ptr, &
                                 rt_sec_transport_vars(ghosted_id), &
                                 reaction,initial_condition, &
                                 sec_tran_constraint,option)
    enddo      
    patch%aux%SC_RT%sec_transport_vars => rt_sec_transport_vars      
  endif

!===============================================================================   

    
  ! allocate aux_var data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(patch%aux%RT%aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RTAuxVarInit(patch%aux%RT%aux_vars(ghosted_id),reaction,option)
  enddo
  patch%aux%RT%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%RT%aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RTAuxVarInit(patch%aux%RT%aux_vars_bc(iconn),reaction,option)
    enddo
  endif
  patch%aux%RT%num_aux_bc = sum_connection

  ! count the number of boundary connections and allocate
  ! aux_var data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%RT%aux_vars_ss(sum_connection))
    do iconn = 1, sum_connection
      call RTAuxVarInit(patch%aux%RT%aux_vars_ss(iconn),reaction,option)
    enddo
  endif
  patch%aux%RT%num_aux_ss = sum_connection
  option%iflag = 0

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RTCreateZeroArray(patch,reaction,option)
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    patch%aux%RT%rt_parameter%diffusion_coefficient(iphase) = &
      cur_fluid_property%diffusion_coefficient
    patch%aux%RT%rt_parameter%diffusion_activation_energy(iphase) = &
      cur_fluid_property%diffusion_activation_energy
    cur_fluid_property => cur_fluid_property%next
  enddo
  
end subroutine RTSetupPatch


! ************************************************************************** !
!
! RTCheckUpdate: In the case of the log formulation, ensures that the update 
!                vector does not exceed a prescribed tolerance
! author: Glenn Hammond
! date: 03/16/09
!
! ************************************************************************** !
#ifndef HAVE_SNES_API_3_2
subroutine RTCheckUpdate(line_search,C,dC,changed,realization,ierr)
#else
subroutine RTCheckUpdate(snes_,C,dC,realization,changed,ierr)
#endif
 
  use Realization_module
 
  implicit none
  
#ifndef HAVE_SNES_API_3_2
  SNESLineSearch :: line_search
#else
  SNES :: snes_
#endif
  Vec :: C
  Vec :: dC
  PetscBool :: changed
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifndef HAVE_SNES_API_3_2
  call RTCheckUpdatePatch(line_search,C,dC,changed,realization,ierr)
#else
  call RTCheckUpdatePatch(snes_,C,dC,realization,changed,ierr)
#endif

end subroutine RTCheckUpdate

! ************************************************************************** !
!
! RTCheckUpdatePatch: In the case of the log formulation, ensures that the 
!                     update vector does not exceed a prescribed tolerance
! author: Glenn Hammond
! date: 03/16/09
!
! ************************************************************************** !
#ifndef HAVE_SNES_API_3_2
subroutine RTCheckUpdatePatch(line_search,C,dC,changed,realization,ierr)
#else
subroutine RTCheckUpdatePatch(snes_,C,dC,realization,changed,ierr)
#endif

  use Realization_module
  use Grid_module
 
  implicit none
  
#ifndef HAVE_SNES_API_3_2
  SNESLineSearch :: line_search
#else
  SNES :: snes_
#endif
  Vec :: C
  Vec :: dC
  type(realization_type) :: realization
  PetscBool :: changed
  
  PetscReal, pointer :: C_p(:)
  PetscReal, pointer :: dC_p(:)
  type(grid_type), pointer :: grid
  PetscReal :: ratio, min_ratio
  PetscInt :: i, n
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  
  call GridVecGetArrayF90(grid,dC,dC_p,ierr)

  if (realization%reaction%use_log_formulation) then
    ! C and dC are actually lnC and dlnC
    dC_p = dsign(1.d0,dC_p)*min(dabs(dC_p),realization%reaction%max_dlnC)

    ! at this point, it does not matter whether "changed" is set to true, since it 
    ! is not check in PETSc.  Thus, I don't want to spend time checking for changes
    ! and performing an allreduce for log formulation.
  
  else
    call VecGetLocalSize(C,n,ierr)
    call GridVecGetArrayF90(grid,C,C_p,ierr)
    
    ! C^p+1 = C^p - dC^p
    ! if dC is positive and abs(dC) larger than C
    ! we need to scale the update
    
    ! compute smallest ratio of C to dC
    min_ratio = 1.d20 ! large number
    do i = 1, n
      if (C_p(i) <= dC_p(i)) then
        ratio = abs(C_p(i)/dC_p(i))
        if (ratio < min_ratio) min_ratio = ratio
      endif
    enddo
    ratio = min_ratio
    
    ! get global minimum
    call MPI_Allreduce(ratio,min_ratio,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                       MPI_MIN,realization%option%mycomm,ierr)
                       
    ! scale if necessary
    if (min_ratio < 1.d0) then
      ! scale by 0.99 to make the update slightly smaller than the min_ratio
      dC_p = dC_p*min_ratio*0.99d0
      changed = PETSC_TRUE
    endif
    call GridVecRestoreArrayF90(grid,C,C_p,ierr)
  endif

  call GridVecRestoreArrayF90(grid,dC,dC_p,ierr)

end subroutine RTCheckUpdatePatch

! ************************************************************************** !
!
! RTComputeMassBalance: 
! author: Glenn Hammond
! date: 12/23/08
!
! ************************************************************************** !
subroutine RTComputeMassBalance(realization,mass_balance)

  use Realization_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%ntrandof, &
                            realization%option%nphase)
  
  mass_balance = 0.d0
  
  call RTComputeMassBalancePatch(realization,mass_balance)

end subroutine RTComputeMassBalance

! ************************************************************************** !
!
! RTComputeMassBalancePatch: Initializes mass balance
! author: Glenn Hammond
! date: 12/23/08
!
! ************************************************************************** !
subroutine RTComputeMassBalancePatch(realization,mass_balance)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%ntrandof, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reaction_type), pointer :: reaction

  PetscReal, pointer :: volume_p(:), porosity_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase
  PetscInt :: i, icomp, imnrl, ncomp, irate, irxn

  iphase = 1
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  reaction => realization%reaction

  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      mass_balance(:,iphase) = mass_balance(:,iphase) + &
        rt_aux_vars(ghosted_id)%total(:,iphase) * &
        global_aux_vars(ghosted_id)%sat(iphase) * &
        porosity_loc_p(ghosted_id) * &
        volume_p(local_id)*1000.d0
        
      if (iphase == 1) then
      ! add contribution of equilibrium sorption
        if (reaction%neqsorb > 0) then
          mass_balance(:,iphase) = mass_balance(:,iphase) + &
            rt_aux_vars(ghosted_id)%total_sorb_eq(:) * volume_p(local_id)
        endif


      ! add contribution of kinetic multirate sorption
        if (reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
          do irxn = 1, reaction%surface_complexation%nkinmrsrfcplxrxn
            do irate = 1, reaction%surface_complexation%kinmr_nrate(irxn)
              mass_balance(:,iphase) = mass_balance(:,iphase) + &
                rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate,irxn) * &
                volume_p(local_id)
            enddo
          enddo
        endif
               
      ! add contribution from mineral volume fractions
        if (reaction%mineral%nkinmnrl > 0) then
          do imnrl = 1, reaction%mineral%nkinmnrl
            ncomp = reaction%mineral%kinmnrlspecid(0,imnrl)
            do i = 1, ncomp
              icomp = reaction%mineral%kinmnrlspecid(i,imnrl)
              mass_balance(icomp,iphase) = mass_balance(icomp,iphase) &
              + reaction%mineral%kinmnrlstoich(i,imnrl)                  &
              * rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl)      &
              * volume_p(local_id) &
              / reaction%mineral%kinmnrl_molar_vol(imnrl)
            enddo 
          enddo
        endif
      endif
    enddo
  enddo

  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine RTComputeMassBalancePatch

! ************************************************************************** !
!
! RTZeroMassBalanceDeltaPatch: Zeros mass balance delta array
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RTZeroMassBalanceDeltaPatch(realization)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  rt_aux_vars_ss => patch%aux%RT%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%RT%num_aux
    patch%aux%RT%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  do iconn = 1, patch%aux%RT%num_aux_bc
    rt_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
  enddo

  do iconn = 1, patch%aux%RT%num_aux_ss
    rt_aux_vars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine RTZeroMassBalanceDeltaPatch

! ************************************************************************** !
!
! RTUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RTUpdateMassBalancePatch(realization)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  rt_aux_vars_ss => patch%aux%RT%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%RT%num_aux
    patch%aux%RT%aux_vars(iconn)%mass_balance = &
      patch%aux%RT%aux_vars(iconn)%mass_balance + &
      patch%aux%RT%aux_vars(iconn)%mass_balance_delta*option%tran_dt
  enddo
#endif

  do iconn = 1, patch%aux%RT%num_aux_bc
    rt_aux_vars_bc(iconn)%mass_balance = &
      rt_aux_vars_bc(iconn)%mass_balance + &
      rt_aux_vars_bc(iconn)%mass_balance_delta*option%tran_dt
  enddo

  do iconn = 1, patch%aux%RT%num_aux_ss
    rt_aux_vars_ss(iconn)%mass_balance = &
      rt_aux_vars_ss(iconn)%mass_balance + &
      rt_aux_vars_ss(iconn)%mass_balance_delta*option%tran_dt
  enddo

end subroutine RTUpdateMassBalancePatch

! ************************************************************************** !
!
! RTInitializeTimestep: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RTInitializeTimestep(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTInitializeTimestepPatch(realization)

end subroutine RTInitializeTimestep

! ************************************************************************** !
!
! RTInitializeTimestepPatch: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTInitializeTimestepPatch(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call RTUpdateFixedAccumulationPatch(realization)
  
  ! geh: never use transport coefs evaluated at time k
!  call RTUpdateTransportCoefsPatch(realization)

end subroutine RTInitializeTimestepPatch
  
! ************************************************************************** !
!
! RTUpdateSolution: Updates data in module after a successful time step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine RTUpdateSolution(realization)

  use Realization_module
  use Discretization_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy,ierr)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%tran_xx, &
                                   realization%field%tran_xx_loc,NTRANDOF)
  
  ! update mineral volume fractions
  call RTUpdateSolutionPatch(realization)

end subroutine RTUpdateSolution

! ************************************************************************** !
!
! RTUpdateSolutionPatch: 
! author: Glenn Hammond
! date: 09/04/08
!
! ************************************************************************** !
subroutine RTUpdateSolutionPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Reaction_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscInt :: ghosted_id, local_id, imnrl, iaqspec, ncomp, icomp
  PetscInt :: k, irate, irxn, icplx, ncplx, ikinrxn
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  PetscReal :: conc, max_conc, min_conc
  PetscErrorCode :: ierr
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars

  ! update:                             cells      bcs         act. coefs.
  call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

!geh: for debugging max/min concentrations
#if 0
  max_conc = -1.d20
  min_conc = 1.d20
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    conc = rt_aux_vars(ghosted_id)%total(1,1)
    max_conc = max(conc,max_conc)
    min_conc = min(conc,min_conc)
  enddo
  call MPI_Allreduce(max_conc,conc,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  max_conc = conc
  call MPI_Allreduce(min_conc,conc,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
  min_conc = conc
  if (option%print_screen_flag) then
    write(*,'("Time: ",1pe12.5," Max: ",1pe12.5," Min: ",1pe12.5)') &
      option%tran_time/realization%output_option%tconv,max_conc, min_conc
  endif
#endif

  if (.not.option%init_stage) then
    ! update mineral volume fractions
    if (reaction%mineral%nkinmnrl > 0) then
    
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle         
        do imnrl = 1, reaction%mineral%nkinmnrl
          ! rate = mol/m^3/sec
          ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
          !                                mol_vol (m^3 mnrl/mol mnrl)
          rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = &
            rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) + &
            rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* &
            reaction%mineral%kinmnrl_molar_vol(imnrl)* &
            option%tran_dt
          if (rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) < 0.d0) &
            rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = 0.d0

#ifdef CHUAN_CO2
          if (option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE) then
            ncomp = reaction%mineral%kinmnrlspecid(0,imnrl)
            do iaqspec=1, ncomp  
              icomp = reaction%mineral%kinmnrlspecid(iaqspec,imnrl)
              if (icomp == realization%reaction%species_idx%co2_aq_id) then
                global_aux_vars(ghosted_id)%reaction_rate(2) &
                  = global_aux_vars(ghosted_id)%reaction_rate(2)& 
                  + rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* option%tran_dt&
                  * reaction%mineral%mnrlstoich(icomp,imnrl)/option%flow_dt
              else if (icomp == reaction%species_idx%h2o_aq_id) then
                global_aux_vars(ghosted_id)%reaction_rate(1) &
                  = global_aux_vars(ghosted_id)%reaction_rate(1)& 
                  + rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* option%tran_dt&
                  * reaction%mineral%mnrlstoich(icomp,imnrl)/option%flow_dt
              endif
            enddo 
          endif   
#endif

        enddo
      enddo
    endif
  

    ! update multirate sorption concentrations 
  ! WARNING: below assumes site concentration multiplicative factor
    if (reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then 
      do ghosted_id = 1, grid%ngmax 
        if (patch%imat(ghosted_id) <= 0) cycle
        do irxn = 1, reaction%surface_complexation%nkinmrsrfcplxrxn
          do irate = 1, reaction%surface_complexation%kinmr_nrate(irxn)
            kdt = reaction%surface_complexation%kinmr_rate(irate,irxn) * &
                  option%tran_dt 
            one_plus_kdt = 1.d0 + kdt 
            k_over_one_plus_kdt = &
              reaction%surface_complexation%kinmr_rate(irate,irxn)/one_plus_kdt
            rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate,irxn) = & 
              (rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate,irxn) + & 
              kdt * reaction%surface_complexation%kinmr_frac(irate,irxn) * &
              rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,0,irxn))/one_plus_kdt
          enddo
        enddo
      enddo 
    endif

    ! update kinetic sorption concentrations
    if (reaction%surface_complexation%nkinsrfcplxrxn > 0) then
      do ghosted_id = 1, grid%ngmax 
        if (patch%imat(ghosted_id) <= 0) cycle
        do ikinrxn = 1, reaction%surface_complexation%nkinsrfcplxrxn
          irxn = reaction%surface_complexation%&
                   kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
          ncplx = reaction%surface_complexation%srfcplxrxn_to_complex(0,irxn)
          do k = 1, ncplx ! ncplx in rxn
            icplx = reaction%surface_complexation%srfcplxrxn_to_complex(k,irxn)
            rt_aux_vars(ghosted_id)%kinsrfcplx_conc(icplx,ikinrxn) = &
              rt_aux_vars(ghosted_id)%kinsrfcplx_conc_kp1(icplx,ikinrxn)
          enddo
        enddo
      enddo
    endif
    
    ! update secondary continuum variables
    if (option%use_mc) then
      do ghosted_id = 1, grid%ngmax
        if (patch%imat(ghosted_id) <= 0) cycle
          sec_diffusion_coefficient = realization% &
                                      material_property_array(1)%ptr% &
                                      secondary_continuum_diff_coeff
          sec_porosity = realization%material_property_array(1)%ptr% &
                         secondary_continuum_porosity

          call SecondaryRTAuxVarComputeMulti(&
                                        rt_sec_transport_vars(ghosted_id), &
                                        rt_aux_vars(ghosted_id), &
                                        global_aux_vars(ghosted_id), &
                                        reaction, &
                                        sec_diffusion_coefficient, &
                                        sec_porosity, &
                                        option)                                        
                                   
      enddo
    endif

  endif

  if (option%compute_mass_balance_new) then
    call RTUpdateMassBalancePatch(realization)
  endif
  
end subroutine RTUpdateSolutionPatch

! ************************************************************************** !
!
! RTUpdateFixedAccumulationPatch: Computes derivative of accumulation term in 
!                           residual function 
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateFixedAccumulationPatch(realization)

  use Realization_module
  use Patch_module
  use Reactive_Transport_Aux_module
  use Option_module
  use Field_module  
  use Grid_module
  use Secondary_Continuum_Aux_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal, pointer :: xx_p(:), porosity_loc_p(:), tor_loc_p(:), &
                        volume_p(:), accum_p(:), density_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: dof_offset, istart, iendaq, iendall
  PetscInt :: istartim, iendim
  PetscInt :: istartcoll, iendcoll
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid
  reaction => realization%reaction
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars

  ! cannot use tran_xx_loc vector here as it has not yet been updated.
  call GridVecGetArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
  
  vol_frac_prim = 1.d0
  
! Do not use RTUpdateAuxVarsPatch() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    ! compute offset in solution vector for first dof in grid cell
    dof_offset = (local_id-1)*reaction%ncomp
    
    ! calculate range of aqueous species
    istart = dof_offset + 1
    iendaq = dof_offset + reaction%naqcomp
    iendall = dof_offset + reaction%ncomp

    ! copy primary aqueous species
    rt_aux_vars(ghosted_id)%pri_molal = xx_p(istart:iendaq)
    
    if (reaction%ncoll > 0) then
      istartcoll = dof_offset + reaction%offset_colloid + 1
      iendcoll = dof_offset + reaction%offset_colloid + reaction%ncoll
      rt_aux_vars(ghosted_id)%colloid%conc_mob = xx_p(istartcoll:iendcoll)* &
        global_aux_vars(ghosted_id)%den_kg(1)*1.d-3
    endif
    
    if (reaction%nimcomp > 0) then
      istartim = dof_offset + reaction%offset_immobile + 1
      iendim = dof_offset + reaction%offset_immobile + reaction%nimcomp
      rt_aux_vars(ghosted_id)%immobile = xx_p(istartim:iendim)
    endif
    
    if (option%use_mc) then
      vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
    endif
    
    ! DO NOT RECOMPUTE THE ACTIVITY COEFFICIENTS BEFORE COMPUTING THE
    ! FIXED PORTION OF THE ACCUMULATION TERM - geh
    call RTAuxVarCompute(rt_aux_vars(ghosted_id), &
                         global_aux_vars(ghosted_id), &
                         reaction,option)
    call RTAccumulation(rt_aux_vars(ghosted_id), &
                        global_aux_vars(ghosted_id), &
                        porosity_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        reaction,option,vol_frac_prim, &
                        accum_p(istart:iendall)) 
    if (reaction%neqsorb > 0) then
      call RAccumulationSorb(rt_aux_vars(ghosted_id), &
                             global_aux_vars(ghosted_id), &
                             volume_p(local_id), &
                             reaction,option,accum_p(istart:iendall))
    endif
  enddo

  call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)

end subroutine RTUpdateFixedAccumulationPatch

! ************************************************************************** !
!
! RTUpdateTransportCoefs: 
! author: Glenn Hammond
! date: 03/09/10
!
! ************************************************************************** !
subroutine RTUpdateTransportCoefs(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTUpdateTransportCoefsPatch(realization)

end subroutine RTUpdateTransportCoefs

! ************************************************************************** !
!
! RTUpdateTransportCoefsPatch: Calculates coefficients for transport matrix 
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
subroutine RTUpdateTransportCoefsPatch(realization)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:)
  PetscInt :: local_id, ghosted_id, ghosted_face_id, id
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set  
  PetscInt :: sum_connection, iconn, num_connections
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter

  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
    ! distance = scalar - magnitude of distance
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TDiffusion(global_aux_vars(ghosted_id_up), &
                      porosity_loc_p(ghosted_id_up), &
                      tor_loc_p(ghosted_id_up), &
                      patch%material_property_array(patch%imat(ghosted_id_up))% &
                        ptr%longitudinal_dispersivity, &
                      dist_up, &
                      global_aux_vars(ghosted_id_dn), &
                      porosity_loc_p(ghosted_id_dn), &
                      tor_loc_p(ghosted_id_dn), &
                      patch%material_property_array(patch%imat(ghosted_id_dn))% &
                        ptr%longitudinal_dispersivity, &
                      dist_dn, &
                      rt_parameter,option, &
                      patch%internal_velocities(:,sum_connection), &
                      patch%internal_tran_coefs(:,sum_connection))
                     
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
 
    if (option%mimetic) then 
      num_connections = boundary_condition%numfaces_set
    else
      cur_connection_set => boundary_condition%connection_set
      num_connections = cur_connection_set%num_connections
    end if
    do iconn = 1, num_connections
      sum_connection = sum_connection + 1
  
      if (option%mimetic) then
#ifdef DASVYAT
        ghosted_face_id = boundary_condition%faces_set(iconn)
        cur_connection_set => grid%faces(ghosted_face_id)%conn_set_ptr
        id = grid%faces(ghosted_face_id)%id
        local_id = grid%nG2L(cur_connection_set%id_dn(id))
#endif
      else
        local_id = cur_connection_set%id_dn(iconn)
      end if
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      call TDiffusionBC(boundary_condition%tran_condition%itype, &
                        global_aux_vars_bc(sum_connection), &
                        global_aux_vars(ghosted_id), &
                        porosity_loc_p(ghosted_id), &
                        tor_loc_p(ghosted_id), &
                        patch%material_property_array(patch%imat(ghosted_id))% &
                          ptr%longitudinal_dispersivity, &
                        cur_connection_set%dist(0,iconn), &
                        rt_parameter,option, &
                        patch%boundary_velocities(:,sum_connection), &
                        patch%boundary_tran_coefs(:,sum_connection))
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)

end subroutine RTUpdateTransportCoefsPatch

! ************************************************************************** !
!
! RTUpdateRHSCoefs: Updates coefficients for the right hand side of linear
!                   transport equation
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTUpdateRHSCoefs(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTUpdateRHSCoefsPatch(realization)

end subroutine RTUpdateRHSCoefs

! ************************************************************************** !
!
! RTUpdateRHSCoefsPatch: Updates coefficients for the right hand side of 
!                        linear transport equation
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTUpdateRHSCoefsPatch(realization)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), rhs_coef_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    rhs_coef_p(local_id) = porosity_loc_p(ghosted_id)* &
                           global_aux_vars(ghosted_id)%sat(iphase)* &
! total already has den_kg within 
!                           global_aux_vars(ghosted_id)%den_kg(iphase)* &
                           1000.d0* &
                           volume_p(local_id)/option%tran_dt
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

end subroutine RTUpdateRHSCoefsPatch

! ************************************************************************** !
!
! RTCalculateRHS_t0: Calculate porition of RHS of transport system
!                  at time t0 or time level k
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateRHS_t0(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTCalculateRHS_t0Patch(realization)

end subroutine RTCalculateRHS_t0

! ************************************************************************** !
!
! RTCalculateRHS_t0Patch: Calculate porition of RHS of transport system
!                         at time t0 or time level k
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateRHS_t0Patch(realization)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: rhs_coef_p(:)
  PetscReal, pointer :: rhs_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscInt :: istartaq, iendaq

  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_aux_vars => patch%aux%RT%aux_vars
  grid => patch%grid
  reaction => realization%reaction

  ! Get vectors
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_rhs,rhs_p,ierr)

  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle    
    iendaq = local_id*reaction%naqcomp
    istartaq = iendaq-reaction%naqcomp+1
    rhs_p(istartaq:iendaq) = rt_aux_vars(ghosted_id)%total(:,iphase)* &
                             rhs_coef_p(local_id)
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tran_rhs,rhs_p,ierr)

end subroutine RTCalculateRHS_t0Patch

! ************************************************************************** !
!
! RTCalculateRHS_t1: Calculate porition of RHS of transport system
!                    at time t1 or time level k+1
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateRHS_t1(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTCalculateRHS_t1Patch(realization)

end subroutine RTCalculateRHS_t1

! ************************************************************************** !
!
! RTCalculateRHS_t1Patch: Calculate porition of RHS of transport system
!                         at time level k+1
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateRHS_t1Patch(realization)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: rhs_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscReal :: coef_up(1), coef_dn(1)
  PetscReal :: msrc(2)
  PetscReal :: Res(realization%reaction%naqcomp)
  PetscInt :: istartaq, iendaq

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: source_sink
  PetscInt :: sum_connection, iconn  
  PetscReal :: qsrc
  PetscInt :: offset, istartcoll, iendcoll, istartall, iendall, icomp, ieqgas
  PetscBool :: volumetric
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in, coef_out
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid
  reaction => realization%reaction

  iphase = 1

!geh - activity coef updates must always be off!!!
!geh    ! update:                             cells      bcs        act. coefs.
!  call RTUpdateAuxVarsPatch(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVarsPatch(realization,PETSC_FALSE,PETSC_TRUE,PETSC_TRUE)
  else
    call RTUpdateAuxVarsPatch(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  endif

  ! Get vectors
  call GridVecGetArrayF90(grid,field%tran_rhs,rhs_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  ! add in inflowing boundary conditions
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

      if (patch%imat(ghosted_id) <= 0) cycle

      call TFluxCoef(option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,sum_connection), &
                     0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                     coef_up,coef_dn)

      ! coef_dn not needed 
      iendaq = local_id*reaction%naqcomp
      istartaq = iendaq-reaction%naqcomp+1
      
      rhs_p(istartaq:iendaq) = rhs_p(istartaq:iendaq) + &
        coef_up(iphase)*rt_aux_vars_bc(sum_connection)%total(:,iphase)

    enddo
    boundary_condition => boundary_condition%next
  enddo  

  ! add in inflowing sources
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    qsrc = 0.d0
    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
    
    ! only handle injection on rhs
    if (qsrc < 0.d0) then
      source_sink => source_sink%next
      cycle
    endif
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_colloid + 1
        iendcoll = reaction%offset_colloid + reaction%ncoll
      endif

      qsrc = patch%ss_fluid_fluxes(1,sum_connection)
      call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)      

      Res(istartaq:iendaq) = & !coef_in*rt_aux_vars(ghosted_id)%total(:,iphase) + &
                             coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                        rt_auxvar%total(:,iphase)
      if (reaction%ncoll > 0) then
        Res(istartcoll:iendcoll) = & !coef_in*rt_aux_vars(ghosted_id)%colloid%conc_mob(:) !+ &
                                   coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                              rt_auxvar%colloid%conc_mob(:)
      endif
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      ! subtract since the contribution is on the rhs
      rhs_p(istartall:iendall) = rhs_p(istartall:iendall) - Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo

#ifdef CHUAN_CO2
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      source_sink => patch%source_sinks%first 
      do 
        if (.not.associated(source_sink)) exit

!geh begin change
!geh        msrc(:) = source_sink%flow_condition%pressure%flow_dataset%time_series%cur_value(:)
        msrc(:) = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(:)
!geh end change
        msrc(1) =  msrc(1) / FMWH2O*1D3
        msrc(2) =  msrc(2) / FMWCO2*1D3
        ! print *,'RT SC source'
        do iconn = 1, cur_connection_set%num_connections      
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          Res=0D0
          
          if (patch%imat(ghosted_id) <= 0) cycle
          
          select case(source_sink%flow_condition%itype(1))
            case(MASS_RATE_SS)
              do ieqgas = 1, reaction%ngas
                if (abs(reaction%species_idx%co2_gas_id) == ieqgas) then
                  icomp = reaction%eqgasspecid(1,ieqgas)
                  iendall = local_id*reaction%ncomp
                  istartall = iendall-reaction%ncomp
                  Res(icomp) = -msrc(2)
                  rhs_p(istartall+icomp) = rhs_p(istartall+icomp) - Res(icomp)
!                 print *,'RT SC source', ieqgas,icomp, res(icomp)  
                endif 
              enddo
          end select 
        enddo
        source_sink => source_sink%next
      enddo
  end select
     
#endif
#endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs,rhs_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

end subroutine RTCalculateRHS_t1Patch

! ************************************************************************** !
!
! RTCalculateTransportMatrix: Calculate transport matrix
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateTransportMatrix(realization,T)

  use Realization_module
  use Option_module
  use Grid_module

  implicit none
      
  type(realization_type) :: realization
  Mat :: T
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option 
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  option => realization%option
 
  call MatZeroEntries(T,ierr)
  
  call RTCalculateTranMatrixPatch1(realization,T)
  call RTCalculateTranMatrixPatch2(realization,T)

  if (realization%debug%matview_Jacobian) then
#if 1
    call PetscViewerASCIIOpen(realization%option%mycomm,'Tmatrix.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Tmatrix.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
#endif
    call MatView(T,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif  
  
end subroutine RTCalculateTransportMatrix

! ************************************************************************** !
!
! RTCalculateTranMatrixPatch: Calculate transport matrix
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateTranMatrixPatch1(realization,T)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  Mat :: T
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: source_sink
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(1), coef_dn(1)
  PetscErrorCode :: ierr
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      call TFluxCoef(option,cur_connection_set%area(iconn), &
                     patch%internal_velocities(:,sum_connection), &
                     patch%internal_tran_coefs(:,sum_connection), &
                     cur_connection_set%dist(-1,iconn), &
                     coef_up,coef_dn)

!      coef_up = coef_up*global_aux_vars(ghosted_id_up)%den_kg*1.d-3
!      coef_dn = coef_dn*global_aux_vars(ghosted_id_dn)%den_kg*1.d-3
                     
      if (local_id_up > 0) then
        call MatSetValuesLocal(T,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr)
        call MatSetValuesLocal(T,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr) 
      endif
      if (local_id_dn > 0) then
        coef_up = -coef_up
        coef_dn = -coef_dn
        call MatSetValuesLocal(T,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr)
        call MatSetValuesLocal(T,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr) 
      endif
                     
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  
  ! add in outflowing boundary conditions
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

      if (patch%imat(ghosted_id) <= 0) cycle

      call TFluxCoef(option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,sum_connection), &
                     0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                     coef_up,coef_dn)

 !     coef_dn = coef_dn*global_aux_vars(ghosted_id)%den_kg*1.d-3

      !Jup not needed 
      coef_dn = -coef_dn
      call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr)
    
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)

  if (patch%aux%RT%inactive_cells_exist) then
    coef = 1.d0
    call MatZeroRowsLocal(T,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,coef, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif
  
end subroutine RTCalculateTranMatrixPatch1


! ************************************************************************** !
!
! RTCalculateTranMatrixPatch: Calculate transport matrix
! author: Glenn Hammond
! date: 04/25/10
!
! ************************************************************************** !
subroutine RTCalculateTranMatrixPatch2(realization,T)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
 
  type(realization_type) :: realization
  Mat :: T
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: source_sink  
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_dn(1)
  PetscReal :: qsrc
  PetscBool :: volumetric  
  PetscErrorCode :: ierr
  PetscInt :: flow_pc
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in, coef_out
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  ! Accumulation term
  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle    
    coef = porosity_loc_p(ghosted_id)* &
           global_aux_vars(ghosted_id)%sat(iphase)* &
!geh           global_aux_vars(ghosted_id)%den_kg(iphase)* &
           1000.d0* &
           volume_p(local_id)/option%tran_dt
    call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef, &
                           ADD_VALUES,ierr)
  enddo
                        
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    qsrc = 0.d0
    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
      
    ! only handle extraction on lhs
    if (qsrc > 0.d0) then
      source_sink => source_sink%next
      cycle
    endif
      
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      qsrc = patch%ss_fluid_fluxes(1,sum_connection)
      call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)

      coef_dn(1) = coef_in
      !geh: do not remove this conditional as otherwise MatSetValuesLocal() 
      !     will be called for injection too (wasted calls)
      if (coef_dn(1) > 0.d0) then
        call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                               ADD_VALUES,ierr)
      endif 

    enddo
    source_sink => source_sink%next
  enddo

  ! All CO2 source/sinks are handled on the RHS for now

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)

  if (patch%aux%RT%inactive_cells_exist) then
    coef = 1.d0
    call MatZeroRowsLocal(T,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,coef, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif
  
end subroutine RTCalculateTranMatrixPatch2

! ************************************************************************** !
!
! RTReact: Calculate reaction
! author: Glenn Hammond
! date: 05/03/10
!
! ************************************************************************** !
subroutine RTReact(realization)

  use Realization_module
  use Field_module
  use Discretization_module    
  use Option_module
  use Logging_module

  implicit none

  type(realization_type) :: realization
  type(discretization_type), pointer :: discretization
  type(field_type), pointer ::field
  type(option_type), pointer :: option

#ifdef OS_STATISTICS
  PetscInt :: call_count
  PetscInt :: sum_newton_iterations
  PetscReal :: ave_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_on_a_core
  PetscInt :: min_newton_iterations_on_a_core
  PetscInt :: temp_int_in(3)
  PetscInt :: temp_int_out(3)
#endif
  
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_rt_react,ierr)
                          
  discretization => realization%discretization
  option => realization%option
  field => realization%field
      
#ifdef OS_STATISTICS
  call_count = 0
  sum_newton_iterations = 0
  max_newton_iterations_in_a_cell = -99999999
  max_newton_iterations_on_a_core = -99999999
  min_newton_iterations_on_a_core = 99999999
#endif  

  call RTReactPatch(realization)

#ifdef OS_STATISTICS
  call_count = call_count + cur_patch%aux%RT%rt_parameter%newton_call_count
  sum_newton_iterations = sum_newton_iterations + &
    cur_patch%aux%RT%rt_parameter%newton_iterations
  if (cur_patch%aux%RT%rt_parameter%max_newton_iterations > &
      max_newton_iterations_in_a_cell) then
    max_newton_iterations_in_a_cell = &
      cur_patch%aux%RT%rt_parameter%max_newton_iterations
  endif
  if (cur_patch%aux%RT%rt_parameter%max_newton_iterations > &
      cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations) then
    cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations = &
      cur_patch%aux%RT%rt_parameter%max_newton_iterations
  endif
#endif 

  ! Logging must come before statistics since the global reductions
  ! will synchonize the cores
  call PetscLogEventEnd(logging%event_rt_react,ierr)
                        
#ifdef OS_STATISTICS
  temp_int_in(1) = call_count
  temp_int_in(2) = sum_newton_iterations
  call MPI_Allreduce(temp_int_in,temp_int_out,TWO_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  ave_newton_iterations_in_a_cell = float(temp_int_out(2)) / temp_int_out(1)

  temp_int_in(1) = max_newton_iterations_in_a_cell
  temp_int_in(2) = sum_newton_iterations ! to calc max # iteration on a core
  temp_int_in(3) = -sum_newton_iterations ! to calc min # iteration on a core
  call MPI_Allreduce(temp_int_in,temp_int_out,THREE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  max_newton_iterations_in_a_cell = temp_int_out(1)
  max_newton_iterations_on_a_core = temp_int_out(2)
  min_newton_iterations_on_a_core = -temp_int_out(3)
  
  if (option%print_screen_flag) then
    write(*, '(" OS Reaction Statistics: ",/, &
             & "   Ave Newton Its / Cell: ",1pe12.4,/, &
             & "   Max Newton Its / Cell: ",i4,/, &
             & "   Max Newton Its / Core: ",i6,/, &
             & "   Min Newton Its / Core: ",i6)') &
               ave_newton_iterations_in_a_cell, &
               max_newton_iterations_in_a_cell, &
               max_newton_iterations_on_a_core, &
               min_newton_iterations_on_a_core
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" OS Reaction Statistics: ",/, &
             & "   Ave Newton Its / Cell: ",1pe12.4,/, &
             & "   Max Newton Its / Cell: ",i4,/, &
             & "   Max Newton Its / Core: ",i6,/, &
             & "   Min Newton Its / Core: ",i6)') &
               ave_newton_iterations_in_a_cell, &
               max_newton_iterations_in_a_cell, &
               max_newton_iterations_on_a_core, &
               min_newton_iterations_on_a_core
  endif

#endif 

end subroutine RTReact

! ************************************************************************** !
!
! RTReactPatch: Calculate reaction
! author: Glenn Hammond
! date: 05/03/10
!
! ************************************************************************** !
subroutine RTReactPatch(realization)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  
  use Secondary_Continuum_Aux_module
  
!$ use omp_lib
     
  implicit none
  
  type(realization_type) :: realization
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend, iendaq
  PetscInt :: iphase
  PetscInt :: ithread, vector_length
  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: mask_p(:)
  PetscInt :: num_iterations
#ifdef OS_STATISTICS
  PetscInt :: sum_iterations
  PetscInt :: max_iterations
#endif
  PetscInt :: icount
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim

  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  rt_aux_vars => patch%aux%RT%aux_vars
  grid => patch%grid
  reaction => realization%reaction
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars

  ! need up update aux vars based on current density/saturation,
  ! but NOT activity coefficients
  call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

  ! Get vectors
  call GridVecGetArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
      
  vol_frac_prim = 1.d0    
  iphase = 1
  ithread = 1
#ifdef OS_STATISTICS
  sum_iterations = 0
  max_iterations = 0
  icount = 0
#endif

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    
    istart = (local_id-1)*reaction%ncomp+1
    iend = istart + reaction%ncomp - 1
    iendaq = istart + reaction%naqcomp - 1
    
    if (option%use_mc) then
      vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
    endif
    
    call RReact(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                tran_xx_p(istart:iend),volume_p(local_id), &
                porosity_loc_p(ghosted_id), &
                num_iterations,reaction,option,vol_frac_prim)
    ! set primary dependent var back to free-ion molality
    tran_xx_p(istart:iendaq) = rt_aux_vars(ghosted_id)%pri_molal
    if (reaction%nimcomp > 0) then
      tran_xx_p(reaction%offset_immobile: &
                reaction%offset_immobile + reaction%nimcomp) = &
        rt_aux_vars(ghosted_id)%immobile
    endif
#ifdef OS_STATISTICS
    if (num_iterations > max_iterations) then
      max_iterations = num_iterations
    endif
    sum_iterations = sum_iterations + num_iterations
    icount = icount + 1
#endif
  enddo
  
#ifdef OS_STATISTICS
  patch%aux%RT%rt_parameter%newton_call_count = icount
  patch%aux%RT%rt_parameter%sum_newton_call_count = &
    patch%aux%RT%rt_parameter%sum_newton_call_count + dble(icount)
  patch%aux%RT%rt_parameter%newton_iterations = sum_iterations
  patch%aux%RT%rt_parameter%sum_newton_iterations = &
    patch%aux%RT%rt_parameter%sum_newton_iterations + dble(sum_iterations)
  patch%aux%RT%rt_parameter%max_newton_iterations = max_iterations
#endif
  
  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  if (option%compute_mass_balance_new) then
    call RTZeroMassBalanceDeltaPatch(realization)
    call RTComputeBCMassBalanceOSPatch(realization)
  endif

end subroutine RTReactPatch

! ************************************************************************** !
!
! RTComputeBCMassBalanceOSPatch: Calculates mass balance at boundary 
!                                conditions for operator split mode
! author: Glenn Hammond
! date: 05/04/10
!
! ************************************************************************** !
subroutine RTComputeBCMassBalanceOSPatch(realization)

  use Realization_module
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none

  type(realization_type) :: realization  

  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:) 
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: flow_src_sink_type
  PetscReal :: qsrc
  
  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)
  PetscReal :: coef_in, coef_out
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  rt_aux_vars_ss => patch%aux%RT%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
  
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

      if (patch%imat(ghosted_id) <= 0) cycle

      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,sum_connection), &
                     0.5d0, &
                     coef_up,coef_dn)
      ! TFlux accomplishes the same as what TBCFlux would
      call TFlux(rt_parameter, &
                 rt_aux_vars_bc(sum_connection), &
                 global_aux_vars_bc(sum_connection), &
                 rt_aux_vars(ghosted_id), &
                 global_aux_vars(ghosted_id), &
                 coef_up,coef_dn,option,Res)

    ! contribution to boundary 
      rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
        rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
!        ! contribution to internal 
!        rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) + Res
    
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set

    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
      
    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      qsrc = patch%ss_fluid_fluxes(1,sum_connection)
      call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)
      
      Res = coef_in*rt_aux_vars(ghosted_id)%total(:,iphase) + &
            coef_out*source_sink%tran_condition%cur_constraint_coupler% &
            rt_auxvar%total(:,iphase)
      if (option%compute_mass_balance_new) then
        ! contribution to boundary 
        rt_aux_vars_ss(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_aux_vars_ss(sum_connection)%mass_balance_delta(:,iphase) + Res
        ! contribution to internal 
!        rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) - Res
        endif
    enddo
    source_sink => source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

end subroutine RTComputeBCMassBalanceOSPatch

! ************************************************************************** !
!
! RTNumericalJacobianTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RTNumericalJacobianTest(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  PetscReal :: derivative, perturbation
  
  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  
  PetscInt :: idof, idof2, icell

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  call VecDuplicate(field%tran_xx,xx_pert,ierr)
  call VecDuplicate(field%tran_xx,res,ierr)
  call VecDuplicate(field%tran_xx,res_pert,ierr)
  
  call MatCreate(option%mycomm,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call RTResidual(PETSC_NULL_OBJECT,field%tran_xx,res,realization,ierr)
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
  do idof = 1,grid%nlmax*option%ntrandof
    icell = (idof-1)/option%ntrandof+1
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    call veccopy(field%tran_xx,xx_pert,ierr)
    call vecgetarrayf90(xx_pert,vec_p,ierr)
    perturbation = vec_p(idof)*perturbation_tolerance
    vec_p(idof) = vec_p(idof)+perturbation
    call vecrestorearrayf90(xx_pert,vec_p,ierr)
    call RTResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
    call vecgetarrayf90(res_pert,vec_p,ierr)
    do idof2 = 1, grid%nlmax*option%ntrandof
      derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
      if (dabs(derivative) > 1.d-30) then
        call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
      endif
    enddo
    call GridVecRestoreArrayF90(grid,res_pert,vec_p,ierr)
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'RTnumerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine RTNumericalJacobianTest

! ************************************************************************** !
!
! RTResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RTResidual(snes,xx,r,realization,ierr)

  use Realization_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscReal, pointer :: xx_p(:), log_xx_p(:)
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscViewer :: viewer  
  
  call PetscLogEventBegin(logging%event_rt_residual,ierr)

  patch => realization%patch
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  if (realization%reaction%use_log_formulation) then
    ! have to convert the log concentration to non-log form
    call GridVecGetArrayF90(patch%grid,field%tran_xx,xx_p,ierr)
    call GridVecGetArrayF90(patch%grid,xx,log_xx_p,ierr)
    xx_p(:) = exp(log_xx_p(:))
    call GridVecRestoreArrayF90(patch%grid,field%tran_xx,xx_p,ierr)
    call GridVecRestoreArrayF90(patch%grid,xx,log_xx_p,ierr)  
    call DiscretizationGlobalToLocal(discretization,field%tran_xx,field%tran_xx_loc,NTRANDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc,NTRANDOF)
  endif

  ! pass #1 for internal and boundary flux terms
  call RTResidualPatch1(snes,xx,r,realization,ierr)

  ! pass #2 for everything else
  call RTResidualPatch2(snes,xx,r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RTresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RTxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call PetscLogEventEnd(logging%event_rt_residual,ierr)

end subroutine RTResidual

! ************************************************************************** !
!
! RTResidualPatch1: Computes residual function for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTResidualPatch1(snes,xx,r,realization,ierr)

  use Realization_module
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Secondary_Continuum_Aux_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim
  
#ifdef CENTRAL_DIFFERENCE  
  PetscReal :: T_11(realization%option%nphase)
  PetscReal :: T_12(realization%option%nphase)
  PetscReal :: T_21(realization%option%nphase)
  PetscReal :: T_22(realization%option%nphase)
  PetscReal :: Res_1(realization%reaction%ncomp)
  PetscReal :: Res_2(realization%reaction%ncomp)
#else
  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)
  PetscReal :: Res(realization%reaction%ncomp)
#endif

#ifdef CHUAN_CO2
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscInt :: icomp, ieqgas
#endif

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars

  
  if (.not.patch%aux%RT%aux_vars_up_to_date) then
    if (reaction%act_coef_update_frequency == ACT_COEF_FREQUENCY_NEWTON_ITER) then
      ! update: cells      bcs        act. coefs.
      call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
    else 
      ! update: cells      bcs        act. coefs.
      call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)
    endif
  endif
  patch%aux%RT%aux_vars_up_to_date = PETSC_FALSE 
  
  if (option%compute_mass_balance_new) then
    call RTZeroMassBalanceDeltaPatch(realization)
  endif
  
  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,r, r_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)

  r_p = 0.d0
  vol_frac_prim = 1.d0

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      ! TFluxCoef will eventually be moved to another routine where it should be
      ! called only once per flux interface at the beginning of a transport
      ! time step.
      
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id_up)%epsilon
      endif  
      
      
#ifndef CENTRAL_DIFFERENCE        
      call TFluxCoef(option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                coef_up,coef_dn)
                      
      call TFlux(rt_parameter, &
                  rt_aux_vars(ghosted_id_up), &
                  global_aux_vars(ghosted_id_up), &
                  rt_aux_vars(ghosted_id_dn), &
                  global_aux_vars(ghosted_id_dn), &
                  coef_up,coef_dn,option,Res)



#ifdef COMPUTE_INTERNAL_MASS_FLUX
      rt_aux_vars(local_id_up)%mass_balance_delta(:,iphase) = &
        rt_aux_vars(local_id_up)%mass_balance_delta(:,iphase) - Res        
#endif
      
      if (local_id_up>0) then
        iend = local_id_up*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)
      endif
      
      if (local_id_dn>0) then
        iend = local_id_dn*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:reaction%ncomp)
      endif

      if (option%store_solute_fluxes) then
        patch%internal_fluxes(iphase,1:reaction%ncomp,iconn) = &
            Res(1:reaction%ncomp)
      endif
#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                 patch%internal_velocities(:,sum_connection), &
                 patch%internal_tran_coefs(:,sum_connection)*vol_frac_prim, &
                 cur_connection_set%dist(-1,iconn), &
                 T_11,T_12,T_21,T_22)
      call TFlux_CD(rt_parameter, &
                  rt_aux_vars(ghosted_id_up), &
                  global_aux_vars(ghosted_id_up), &
                  rt_aux_vars(ghosted_id_dn), &
                  global_aux_vars(ghosted_id_dn), &
                  T_11,T_12,T_21,T_22,option,Res_1,Res_2)
                             
      if (local_id_up>0) then
        iend = local_id_up*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res_1(1:reaction%ncomp)
      endif
      
      if (local_id_dn>0) then
        iend = local_id_dn*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res_2(1:reaction%ncomp)
      endif
#endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo
    
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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
      endif  
      
#ifndef CENTRAL_DIFFERENCE
      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(option,cur_connection_set%area(iconn), &
                  patch%boundary_velocities(:,sum_connection), &
                  patch%boundary_tran_coefs(:,sum_connection)*vol_frac_prim, &
                  0.5d0, &
                  coef_up,coef_dn)
      ! TFlux accomplishes the same as what TBCFlux would
      call TFlux(rt_parameter, &
                  rt_aux_vars_bc(sum_connection), &
                  global_aux_vars_bc(sum_connection), &
                  rt_aux_vars(ghosted_id), &
                  global_aux_vars(ghosted_id), &
                  coef_up,coef_dn,option,Res)
                  
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:reaction%ncomp)

      if (option%store_solute_fluxes) then
        patch%boundary_fluxes(iphase,1:reaction%ncomp,sum_connection) = &
            -Res(1:reaction%ncomp)
      endif

      if (option%compute_mass_balance_new) then
      ! contribution to boundary 
        rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
!        ! contribution to internal 
!        rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) + Res
        endif  

#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                patch%boundary_velocities(:,sum_connection), &
                patch%boundary_tran_coefs(:,sum_connection)*vol_frac_prim, &
                0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                T_11,T_12,T_21,T_22)
      call TFlux_CD(rt_parameter, &
                  rt_aux_vars_bc(sum_connection), &
                  global_aux_vars_bc(sum_connection), &
                  rt_aux_vars(ghosted_id), &
                  global_aux_vars(ghosted_id), &
                  T_11,T_12,T_21,T_22,option,Res_1,Res_2)

      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) + Res_2(1:reaction%ncomp)

      if (option%store_solute_fluxes) then
        patch%boundary_fluxes(iphase,1:reaction%ncomp,sum_connection) = &
            Res_2(1:reaction%ncomp)
      endif
      if (option%compute_mass_balance_new) then
      ! contribution to boundary 
        rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_aux_vars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res_2
!        ! contribution to internal 
!        rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) + Res
        endif  
      
#endif                   
     
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)

end subroutine RTResidualPatch1

! ************************************************************************** !
!
! RTResidualPatch2: Computes residual function for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTResidualPatch2(snes,xx,r,realization,ierr)

  use Realization_module
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Logging_module
  use Secondary_Continuum_Aux_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  PetscInt :: i
  PetscInt :: istartaq, iendaq
  PetscInt :: istartcoll, iendcoll
  PetscInt :: istartall, iendall
  PetscInt :: idof
  PetscInt :: offset
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscReal :: qsrc, molality
  PetscInt :: flow_src_sink_type
  PetscReal :: scale, coef_in, coef_out
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscBool :: volumetric
  PetscInt :: sum_connection
#ifdef CHUAN_CO2
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscInt :: icomp, ieqgas
#endif

  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: res_sec_transport(realization%reaction%ncomp)

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_ss => patch%aux%RT%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  
  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  
  vol_frac_prim = 1.d0

  if (.not.option%steady_state) then
    r_p = r_p - accum_p
#if 1
    ! Accumulation terms ------------------------------------
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
      endif  
        
      call RTAccumulation(rt_aux_vars(ghosted_id), &
                          global_aux_vars(ghosted_id), &
                          porosity_loc_p(ghosted_id), &
                          volume_p(local_id), &
                          reaction,option,vol_frac_prim,Res)
      if (reaction%neqsorb > 0) then
        call RAccumulationSorb(rt_aux_vars(ghosted_id), &
                               global_aux_vars(ghosted_id), &
                               volume_p(local_id), &
                               reaction,option,Res)
      endif
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)
      if (reaction%calculate_water_age) then 
        call RAge(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                  porosity_loc_p(ghosted_id),volume_p(local_id),option,reaction,Res)
        r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)
      endif
      if (reaction%calculate_tracer_age) then 
        call RAge(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                  porosity_loc_p(ghosted_id),volume_p(local_id),option,reaction,Res)
        r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)
      endif
    enddo
  endif
#endif
#if 1

! ========== Secondary continuum transport source terms -- MULTICOMPONENT ======
  if (option%use_mc) then
  ! Secondary continuum contribution (SK 1/31/2012)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp
         
      sec_diffusion_coefficient = realization% &
                                  material_property_array(1)%ptr% &
                                  secondary_continuum_diff_coeff
      sec_porosity = realization%material_property_array(1)%ptr% &
                     secondary_continuum_porosity

      call RTSecondaryTransportMulti(rt_sec_transport_vars(ghosted_id), &
                                     rt_aux_vars(ghosted_id), &
                                     global_aux_vars(ghosted_id), &
                                     reaction, &
                                     sec_diffusion_coefficient, &
                                     sec_porosity, &
                                     option,res_sec_transport)

      r_p(istartall:iendall) = r_p(istartall:iendall) - &
                               res_sec_transport(1:reaction%ncomp)&
                               *volume_p(local_id)*1.d3 ! convert vol to L from m3

    enddo   
  endif
! ============== end secondary continuum coupling terms ========================

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_colloid + 1
        iendcoll = reaction%offset_colloid + reaction%ncoll
      endif

      qsrc = patch%ss_fluid_fluxes(1,sum_connection)
      call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)

      
      Res(istartaq:iendaq) = coef_in*rt_aux_vars(ghosted_id)%total(:,iphase) + &
                             coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                        rt_auxvar%total(:,iphase)
      
      if (reaction%ncoll > 0) then
        Res(istartcoll:iendcoll) = coef_in*rt_aux_vars(ghosted_id)%colloid%conc_mob(:) + &
                                   coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                              rt_auxvar%colloid%conc_mob(:)
      endif
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)                                  
      if (option%compute_mass_balance_new) then
        ! contribution to boundary 
        rt_aux_vars_ss(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_aux_vars_ss(sum_connection)%mass_balance_delta(:,iphase) + Res
        ! contribution to internal 
!        rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars(ghosted_id)%mass_balance_delta(:,iphase) - Res
        endif
    enddo
    source_sink => source_sink%next
  enddo

#ifdef CHUAN_CO2
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      source_sink => patch%source_sinks%first 
      do 
        if (.not.associated(source_sink)) exit

        select case(source_sink%flow_condition%itype(1))
          case(MASS_RATE_SS)
            msrc(:) = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(:)
          case default
            msrc(:) = 0.d0
        end select

        msrc(1) =  msrc(1) / FMWH2O*1D3
        msrc(2) =  msrc(2) / FMWCO2*1D3
        ! print *,'RT SC source'
        do iconn = 1, cur_connection_set%num_connections      
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          Res=0D0
          
          if (patch%imat(ghosted_id) <= 0) cycle
          
          select case(source_sink%flow_condition%itype(1))
            case(MASS_RATE_SS)
              do ieqgas = 1, reaction%ngas
                if (abs(reaction%species_idx%co2_gas_id) == ieqgas) then
                  icomp = reaction%eqgasspecid(1,ieqgas)
                  iendall = local_id*reaction%ncomp
                  istartall = iendall-reaction%ncomp
                  Res(icomp) = -msrc(2)
                  r_p(istartall+icomp) = r_p(istartall+icomp) + Res(icomp)
!                 print *,'RT SC source', ieqgas,icomp, res(icomp)  
                endif 
              enddo
          end select 
        enddo
        source_sink => source_sink%next
      enddo
  end select
     
#endif
#endif
#if 1  
! Reactions
  if (associated(reaction)) then
  
    call PetscLogEventBegin(logging%event_rt_res_reaction,ierr)  
    
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      Res = 0.d0
      Jup = 0.d0
      call RReaction(Res,Jup,PETSC_FALSE,rt_aux_vars(ghosted_id), &
                     global_aux_vars(ghosted_id), &
                     porosity_loc_p(ghosted_id), &
                     volume_p(local_id),reaction,option)
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
        Res = Res*vol_frac_prim
      endif 
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)                    

    enddo

    call PetscLogEventEnd(logging%event_rt_res_reaction,ierr)   
  endif
#endif

  if (patch%aux%RT%inactive_cells_exist) then
    do i=1,patch%aux%RT%n_zero_rows
      r_p(patch%aux%RT%zero_rows_local(i)) = 0.d0
    enddo
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

end subroutine RTResidualPatch2

! ************************************************************************** !
!
! RTJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RTJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Logging_module

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
  type(grid_type),  pointer :: grid

  call PetscLogEventBegin(logging%event_rt_jacobian,ierr)

#if 0
  call RTNumericalJacobianTest(realization)
#endif

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


  call PetscLogEventBegin(logging%event_rt_jacobian1,ierr)


  ! pass #1 for internal and boundary flux terms  
  call RTJacobianPatch1(snes,xx,J,J,flag,realization,ierr)

  call PetscLogEventEnd(logging%event_rt_jacobian1,ierr)
  call PetscLogEventBegin(logging%event_rt_jacobian2,ierr)
  
  ! pass #2 for everything else
  call RTJacobianPatch2(snes,xx,J,J,flag,realization,ierr)

  call PetscLogEventEnd(logging%event_rt_jacobian2,ierr)
    
  if (realization%debug%matview_Jacobian) then
#if 1
    call PetscViewerASCIIOpen(realization%option%mycomm,'RTjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'RTjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
#endif
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  if (realization%reaction%use_log_formulation) then
    call MatDiagonalScaleLocal(J,realization%field%tran_work_loc,ierr)

    if (realization%debug%matview_Jacobian) then
#if 1
      call PetscViewerASCIIOpen(realization%option%mycomm,'RTjacobianLog.out', &
                                viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'RTjacobianLog.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
#endif
      call MatView(J,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)
    endif
    
  endif

  call PetscLogEventEnd(logging%event_rt_jacobian,ierr)
  
end subroutine RTJacobian

! ************************************************************************** !
!
! RTJacobianPatch1: Computes Jacobian for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTJacobianPatch1(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Logging_module  
  use Secondary_Continuum_Aux_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag  
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reaction_type), pointer :: reaction
      
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim

#ifdef CENTRAL_DIFFERENCE
  PetscReal :: T_11(realization%option%nphase)
  PetscReal :: T_12(realization%option%nphase)
  PetscReal :: T_21(realization%option%nphase)
  PetscReal :: T_22(realization%option%nphase)
  PetscReal :: J_11(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_12(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_21(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_22(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
#else
  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Jdn(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
#endif

  option => realization%option
  field => realization%field
  patch => realization%patch  
  grid => patch%grid
  reaction => realization%reaction
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars


  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)
  
  vol_frac_prim = 1.d0

  ! Interior Flux Terms -----------------------------------
  ! must zero out Jacobian blocks

  call PetscLogEventBegin(logging%event_rt_jacobian_flux,ierr)

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id_up)%epsilon
      endif 

#ifndef CENTRAL_DIFFERENCE
      call TFluxCoef(option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                coef_up,coef_dn)
      call TFluxDerivative(rt_parameter, &
                           rt_aux_vars(ghosted_id_up), &
                           global_aux_vars(ghosted_id_up), &
                           rt_aux_vars(ghosted_id_dn), &
                           global_aux_vars(ghosted_id_dn), &
                           coef_up,coef_dn,option,Jup,Jdn)
      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)        
      endif
   
      if (local_id_dn>0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
      endif

#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                T_11,T_12,T_21,T_22)
      call TFluxDerivative_CD(rt_parameter, &
                           rt_aux_vars(ghosted_id_up), &
                           global_aux_vars(ghosted_id_up), &
                           rt_aux_vars(ghosted_id_dn), &
                           global_aux_vars(ghosted_id_dn), &
                           T_11,T_12,T_21,T_22,option, &
                           J_11,J_12,J_21,J_22)
      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      J_11,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      J_12,ADD_VALUES,ierr)        
      endif
   
      if (local_id_dn>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      J_22,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      J_21,ADD_VALUES,ierr)
      endif
#endif


    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  call PetscLogEventEnd(logging%event_rt_jacobian_flux,ierr)
  
  ! Boundary Flux Terms -----------------------------------
  ! must zero out Jacobian block

  call PetscLogEventBegin(logging%event_rt_jacobian_fluxbc,ierr)

  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
    
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
      endif 

#ifndef CENTRAL_DIFFERENCE
      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(option,cur_connection_set%area(iconn), &
                patch%boundary_velocities(:,sum_connection), &
                patch%boundary_tran_coefs(:,sum_connection)*vol_frac_prim, &
                0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                coef_up,coef_dn)
      ! TFluxDerivative accomplishes the same as what TBCFluxDerivative would
      call TFluxDerivative(rt_parameter, &
                           rt_aux_vars_bc(sum_connection), &
                           global_aux_vars_bc(sum_connection), &
                           rt_aux_vars(ghosted_id), &
                           global_aux_vars(ghosted_id), &
                           coef_up,coef_dn,option,Jup,Jdn)

      !Jup not needed 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                 patch%boundary_velocities(:,sum_connection), &
                 patch%boundary_tran_coefs(:,sum_connection)*vol_frac_prim, &
                 0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                 T_11,T_12,T_21,T_22)
      call TFluxDerivative_CD(rt_parameter, &
                           rt_aux_vars_bc(sum_connection), &
                           global_aux_vars_bc(sum_connection), &
                           rt_aux_vars(ghosted_id), &
                           global_aux_vars(ghosted_id), &
                           T_11,T_12,T_21,T_22,option, &
                           J_11,J_12,J_21,J_22)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,J_22,ADD_VALUES,ierr)
#endif
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  call PetscLogEventEnd(logging%event_rt_jacobian_fluxbc,ierr)  

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)

end subroutine RTJacobianPatch1

! ************************************************************************** !
!
! RTJacobianPatch2: Computes Jacobian for reactive transport
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RTJacobianPatch2(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Logging_module
  use Secondary_Continuum_Aux_module

  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag  
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), work_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istartaq, iendaq
  PetscInt :: istartcoll, iendcoll
  PetscInt :: offset, idof                  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: tran_pc
    
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)    
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: qsrc, rdum
  PetscBool :: volumetric
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in, coef_out
  PetscReal :: scale
  
  ! secondary continuum variables
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: jac_transport(realization%reaction%naqcomp,realization%reaction%naqcomp)
  PetscInt :: ncomp

  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars

  vol_frac_prim = 1.d0
  
  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
    
  if (.not.option%steady_state) then
  call PetscLogEventBegin(logging%event_rt_jacobian_accum,ierr)  
#if 1  
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      if (option%use_mc) then    
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
      endif      
      
      call RTAccumulationDerivative(rt_aux_vars(ghosted_id), &
                                    global_aux_vars(ghosted_id), &
                                    porosity_loc_p(ghosted_id), &
                                    volume_p(local_id),reaction,option, &
                                    vol_frac_prim,Jup) 
                                    
      if (reaction%neqsorb > 0) then
        call RAccumulationSorbDerivative(rt_aux_vars(ghosted_id), &
                                         global_aux_vars(ghosted_id), &
                                         volume_p(local_id),reaction, &
                                         option,Jup)
      endif
      
      if (option%use_mc) then

        sec_diffusion_coefficient = realization%material_property_array(1)% &
                                    ptr%secondary_continuum_diff_coeff
        sec_porosity = realization%material_property_array(1)%ptr% &
                       secondary_continuum_porosity
                        
        if (realization%reaction%ncomp /= realization%reaction%naqcomp) then
          option%io_buffer = 'Current multicomponent implementation is for '// &
                             'aqueous reactions only'
          call printErrMsg(option)
        endif
        
        if (rt_sec_transport_vars(ghosted_id)%sec_jac_update) then
          jac_transport = rt_sec_transport_vars(ghosted_id)%sec_jac
        else
          option%io_buffer = 'RT secondary continuum term in primary '// &
                             'jacobian not updated'
          call printErrMsg(option)
        endif
         
        Jup = Jup - jac_transport*volume_p(local_id)*1.d3     ! convert m3 to L                                                                    
                                                                                
      endif

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
    enddo
#endif
  call PetscLogEventEnd(logging%event_rt_jacobian_accum,ierr)  
  endif
#if 1
  ! Source/Sink terms -------------------------------------
  call PetscLogEventBegin(logging%event_rt_jacobian_ss,ierr)   
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

      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_colloid + 1
        iendcoll = reaction%offset_colloid + reaction%ncoll
      endif

      qsrc = patch%ss_fluid_fluxes(1,sum_connection)
      call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype,coef_in,coef_out)

      ! coef_in is non-zero
      if (dabs(coef_in-1.d20) > 0.d0) then
        Jup = coef_in*rt_aux_vars(ghosted_id)%aqueous%dtotal(:,:,option%liquid_phase)         
        if (reaction%ncoll > 0) then
          option%io_buffer = 'Source/sink not yet implemented for colloids'
          call printErrMsg(option)
        endif
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr) 
      endif
    enddo                       
    source_sink => source_sink%next
  enddo
  
  call PetscLogEventEnd(logging%event_rt_jacobian_ss,ierr)  
#endif


#if 1  
! Reactions
  if (associated(reaction)) then

    call PetscLogEventBegin(logging%event_rt_jac_reaction,ierr)
                              
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      Res = 0.d0
      Jup = 0.d0
      call RReactionDerivative(Res,Jup,rt_aux_vars(ghosted_id), &
                               global_aux_vars(ghosted_id), &
                               porosity_loc_p(ghosted_id), &
                               volume_p(local_id),reaction,option)
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
        Jup = Jup*vol_frac_prim
      endif
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1, &
                                    Jup,ADD_VALUES,ierr)                        
    enddo
    
    call PetscLogEventEnd(logging%event_rt_jac_reaction,ierr)
    
  endif
#endif
 
  if (reaction%use_log_formulation) then
    call PetscLogEventBegin(logging%event_rt_jacobian_zero_calc,ierr)  
    call GridVecGetArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      offset = (ghosted_id-1)*reaction%ncomp
      istartaq = offset + reaction%offset_aqueous + 1
      iendaq = offset + reaction%offset_aqueous + reaction%naqcomp
      if (reaction%ncoll > 0) then
        istartcoll = offset + reaction%offset_colloid + 1
        iendcoll = offset + reaction%offset_colloid + reaction%ncoll
      endif
      if (patch%imat(ghosted_id) <= 0) then
        work_loc_p(istartaq:iendaq) = 1.d0
        if (reaction%ncoll > 0) then
          work_loc_p(istartcoll:iendcoll) = 1.d0
        endif
      else
        work_loc_p(istartaq:iendaq) = rt_aux_vars(ghosted_id)%pri_molal(:)
        if (reaction%ncoll > 0) then
          work_loc_p(istartcoll:iendcoll) = rt_aux_vars(ghosted_id)%colloid%conc_mob(:)
        endif
      endif
    enddo
    call GridVecRestoreArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
    call PetscLogEventEnd(logging%event_rt_jacobian_zero_calc,ierr)    
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  
  if (patch%aux%RT%inactive_cells_exist) then
    call PetscLogEventBegin(logging%event_rt_jacobian_zero,ierr)    
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,rdum, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
    call PetscLogEventEnd(logging%event_rt_jacobian_zero,ierr)                          
  endif

end subroutine RTJacobianPatch2

! ************************************************************************** !
!
! RTUpdateAuxVars: Updates the auxiliary variables associated with 
!                  reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVars(realization,update_cells,update_bcs, &
                           update_activity_coefs)

  use Realization_module

  type(realization_type) :: realization
  PetscBool :: update_cells
  PetscBool :: update_bcs
  PetscBool :: update_activity_coefs
  
  ! PETSC_TRUE to update cells
  call RTUpdateAuxVarsPatch(realization,update_cells,update_bcs, &
                            update_activity_coefs)
  
end subroutine RTUpdateAuxVars

! ************************************************************************** !
!
! RTUpdateAuxVarsPatch: Updates the auxiliary variables associated with 
!                       reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVarsPatch(realization,update_cells,update_bcs, &
                                compute_activity_coefs)

  use Realization_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  use Logging_module
  
  implicit none

  type(realization_type) :: realization
  PetscBool :: update_bcs
  PetscBool :: update_cells
  PetscBool :: compute_activity_coefs
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: istartaq, iendaq 
  PetscInt :: istartcoll, iendcoll
  PetscInt :: istartaq_loc, iendaq_loc
  PetscInt :: istartcoll_loc, iendcoll_loc
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal :: xxbc(realization%reaction%ncomp)
  PetscReal, pointer :: basis_molarity_p(:)
  PetscReal, pointer :: basis_coll_conc_p(:)
  PetscReal :: weight
  PetscInt, parameter :: iphase = 1
  PetscInt :: offset
  PetscErrorCode :: ierr
  PetscBool :: skip_equilibrate_constraint
  PetscInt, save :: icall
  
  data icall/0/

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  reaction => realization%reaction
  
  call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  if (update_cells) then

    call PetscLogEventBegin(logging%event_rt_auxvars,ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials

      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*reaction%ncomp
      istartaq = offset + reaction%offset_aqueous + 1
      iendaq = offset + reaction%offset_aqueous + reaction%naqcomp
      
      patch%aux%RT%aux_vars(ghosted_id)%pri_molal = xx_loc_p(istartaq:iendaq)
      if (reaction%ncoll > 0) then
        istartcoll = offset + reaction%offset_colloid + 1
        iendcoll = offset + reaction%offset_colloid + reaction%ncoll
        patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_mob = xx_loc_p(istartcoll:iendcoll)* &
          patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)*1.d-3
      endif
      
      if (compute_activity_coefs) then
        call RActivityCoefficients(patch%aux%RT%aux_vars(ghosted_id), &
                                   patch%aux%Global%aux_vars(ghosted_id), &
                                   reaction,option)
        if (option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE) then
          call CO2AqActCoeff(patch%aux%RT%aux_vars(ghosted_id), &
                                   patch%aux%Global%aux_vars(ghosted_id), &
                                   reaction,option)
        endif                           
      endif
      call RTAuxVarCompute(patch%aux%RT%aux_vars(ghosted_id), &
                           patch%aux%Global%aux_vars(ghosted_id), &
                           reaction,option)
      if (associated(reaction%species_idx) .and. &
          associated(patch%aux%Global%aux_vars(ghosted_id)%m_nacl)) then
        if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
          patch%aux%Global%aux_vars(ghosted_id)%m_nacl(1) = &
                patch%aux%RT%aux_vars(ghosted_id)%pri_molal(reaction%species_idx%na_ion_id)
          patch%aux%Global%aux_vars(ghosted_id)%m_nacl(2) = &
                patch%aux%RT%aux_vars(ghosted_id)%pri_molal(reaction%species_idx%cl_ion_id)
         else
          patch%aux%Global%aux_vars(ghosted_id)%m_nacl = option%m_nacl
        endif
      endif
    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars,ierr)
  endif

  if (update_bcs) then

    call PetscLogEventBegin(logging%event_rt_auxvars_bc,ierr)

    boundary_condition => patch%boundary_conditions%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      basis_molarity_p => boundary_condition%tran_condition% &
        cur_constraint_coupler%aqueous_species%basis_molarity
        

      if (reaction%ncoll > 0) then
        basis_coll_conc_p => boundary_condition%tran_condition% &
                             cur_constraint_coupler%colloids%basis_conc_mob
      endif

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        
        if (patch%imat(ghosted_id) <= 0) cycle

        offset = (ghosted_id-1)*reaction%ncomp
        istartaq_loc = reaction%offset_aqueous + 1
        iendaq_loc = reaction%offset_aqueous + reaction%naqcomp
        istartaq = offset + istartaq_loc
        iendaq = offset + iendaq_loc
    
        if (reaction%ncoll > 0) then
          istartcoll_loc = reaction%offset_colloid + 1
          iendcoll_loc = reaction%offset_colloid + reaction%ncoll
          istartcoll = offset + istartcoll_loc
          iendcoll = offset + iendcoll_loc
        endif

!       if (option%iflowmode /= MPH_MODE .or. icall>1) then
        if (option%iflowmode /= MPH_MODE .and. option%iflowmode /= FLASH2_MODE) then
  !       Note: the  DIRICHLET_BC is not time dependent in this case (icall)    
          select case(boundary_condition%tran_condition%itype)
            case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
              ! since basis_molarity is in molarity, must convert to molality
                ! by dividing by density of water (mol/L -> mol/kg)
              xxbc(istartaq_loc:iendaq_loc) = &
                basis_molarity_p(1:reaction%naqcomp) / &
                patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                1000.d0
              if (reaction%ncoll > 0) then
                xxbc(istartcoll_loc:iendcoll_loc) = &
                  basis_coll_conc_p(1:reaction%ncoll) / &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                  1000.d0
              endif
            case(DIRICHLET_ZERO_GRADIENT_BC)
  !geh            do iphase = 1, option%nphase
                if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! same as dirichlet above
                  xxbc(istartaq_loc:iendaq_loc) = &
                    basis_molarity_p(1:reaction%naqcomp) / &
                    patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                  & 1000.d0
                  if (reaction%ncoll > 0) then
                    xxbc(istartcoll_loc:iendcoll_loc) = &
                      basis_coll_conc_p(1:reaction%ncoll) / &
                      patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                      1000.d0
                  endif
                else
                  ! same as zero_gradient below
                  xxbc(istartaq_loc:iendaq_loc) = xx_loc_p(istartaq:iendaq)
                  if (reaction%ncoll > 0) then
                    xxbc(istartcoll_loc:iendcoll_loc) = &
                      basis_coll_conc_p(1:reaction%ncoll) / &
                      patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                      1000.d0
                  endif
                endif
  !geh          enddo
            case(ZERO_GRADIENT_BC)
              xxbc(istartaq_loc:iendaq_loc) = xx_loc_p(istartaq:iendaq)
              if (reaction%ncoll > 0) then
                xxbc(istartcoll_loc:iendcoll_loc) = &
                  basis_coll_conc_p(1:reaction%ncoll) / &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * &
                  1000.d0
              endif
          end select
          ! no need to update boundary fluid density since it is already set
          patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = &
            xxbc(istartaq_loc:iendaq_loc)
          if (reaction%ncoll > 0) then
            patch%aux%RT%aux_vars_bc(sum_connection)%colloid%conc_mob = &
              xxbc(istartcoll_loc:iendcoll_loc)* &
              patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(1)*1.d-3
          endif
          if (compute_activity_coefs) then
            call RActivityCoefficients(patch%aux%RT%aux_vars_bc(sum_connection), &
                                        patch%aux%Global%aux_vars_bc(sum_connection), &
                                        reaction,option)
            if (option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE) then
              call CO2AqActCoeff(patch%aux%RT%aux_vars_bc(sum_connection), &
                                  patch%aux%Global%aux_vars_bc(sum_connection), &
                                  reaction,option) 
              endif                           
          endif
          call RTAuxVarCompute(patch%aux%RT%aux_vars_bc(sum_connection), &
                                patch%aux%Global%aux_vars_bc(sum_connection), &
                                reaction,option)
        else
          skip_equilibrate_constraint = PETSC_FALSE
        ! Chuan needs to fill this in.
          select case(boundary_condition%tran_condition%itype)
            case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
              ! don't need to do anything as the constraint below provides all
              ! the concentrations, etc.
              
              !geh: terrible kludge, but should work for now.
              !geh: the problem is that ...%pri_molal() on first call is zero and 
              !     PETSC_TRUE is passed into ReactionEquilibrateConstraint() below
              !     for use_prev_soln_as_guess.  If the previous solution is zero,
              !     the code will crash.
              if (patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal(1) < 1.d-200) then
                patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = 1.d-9
              endif
            case(DIRICHLET_ZERO_GRADIENT_BC)
                if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! don't need to do anything as the constraint below provides all
                  ! the concentrations, etc.
                  
                if (patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal(1) < 1.d-200) then
                  patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = 1.d-9
                endif                  
                else
                  ! same as zero_gradient below
                  skip_equilibrate_constraint = PETSC_TRUE
                  patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = &
                    xx_loc_p(istartaq:iendaq)
                  if (reaction%ncoll > 0) then
                    patch%aux%RT%aux_vars_bc(sum_connection)%colloid%conc_mob = &
                      xx_loc_p(istartcoll:iendcoll)* &
                      patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(1)*1.d-3
                  endif                  
                endif
            case(ZERO_GRADIENT_BC)
              skip_equilibrate_constraint = PETSC_TRUE
              patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = &
                xx_loc_p(istartaq:iendaq)
              if (reaction%ncoll > 0) then
                patch%aux%RT%aux_vars_bc(sum_connection)%colloid%conc_mob = &
                  xx_loc_p(istartcoll:iendcoll)* &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(1)*1.d-3
              endif                
          end select
          ! no need to update boundary fluid density since it is already set
          if (.not.skip_equilibrate_constraint) then
            ! print *,'RT redo constrain on BCs: 1: ', sum_connection
            call ReactionEquilibrateConstraint(patch%aux%RT%aux_vars_bc(sum_connection), &
              patch%aux%Global%aux_vars_bc(sum_connection),reaction, &
              boundary_condition%tran_condition%cur_constraint_coupler%constraint_name, &
              boundary_condition%tran_condition%cur_constraint_coupler%aqueous_species, &
              boundary_condition%tran_condition%cur_constraint_coupler%minerals, &
              boundary_condition%tran_condition%cur_constraint_coupler%surface_complexes, &
              boundary_condition%tran_condition%cur_constraint_coupler%colloids, &
              boundary_condition%tran_condition%cur_constraint_coupler%biomass, &
              porosity_loc_p(ghosted_id), &
              boundary_condition%tran_condition%cur_constraint_coupler%num_iterations, &
              PETSC_TRUE,option)
            ! print *,'RT redo constrain on BCs: 2: ', sum_connection  
          endif         
        endif

        if (associated(reaction%species_idx) .and. &
            associated(patch%aux%Global%aux_vars_bc(sum_connection)%m_nacl)) then
          if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
            patch%aux%Global%aux_vars_bc(sum_connection)%m_nacl(1) = &
                  patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal(reaction%species_idx%na_ion_id)
            patch%aux%Global%aux_vars_bc(sum_connection)%m_nacl(2) = &
                  patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal(reaction%species_idx%cl_ion_id)
            else
            patch%aux%Global%aux_vars_bc(sum_connection)%m_nacl = option%m_nacl
          endif
        endif
      enddo ! iconn
      boundary_condition => boundary_condition%next
    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars_bc,ierr)

  endif 

  patch%aux%RT%aux_vars_up_to_date = update_cells .and. update_bcs
  
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  icall = icall+ 1

end subroutine RTUpdateAuxVarsPatch

! ************************************************************************** !
!
! RTCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RTCreateZeroArray(patch,reaction,option)

  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

  type(patch_type) :: patch
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id, icomp

  type(grid_type), pointer :: grid
  PetscInt :: flag
  PetscInt :: ndof
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr

  flag = 0
  grid => patch%grid
  
  n_zero_rows = 0
  
  if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
    ndof = reaction%ncomp
  else
    ndof = 1
  endif

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      n_zero_rows = n_zero_rows + ndof
    else
    endif
  enddo

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))

  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      do icomp = 1, ndof
        ncount = ncount + 1
        zero_rows_local(ncount) = (local_id-1)*ndof+icomp
        zero_rows_local_ghosted(ncount) = (ghosted_id-1)*ndof+icomp-1
      enddo
    else
    endif
  enddo

  patch%aux%RT%zero_rows_local => zero_rows_local
  patch%aux%RT%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%RT%n_zero_rows = n_zero_rows  

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
     
  if (flag > 0) patch%aux%RT%inactive_cells_exist = PETSC_TRUE
     
  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
   endif

end subroutine RTCreateZeroArray

! ************************************************************************** !
!
! RTMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTMaxChange(realization)

  use Realization_module
  use Option_module
  use Field_module
  use Patch_module
  use Grid_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  PetscReal, pointer :: dxx_ptr(:), xx_ptr(:), yy_ptr(:)
  
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field

  option%dcmax=0.D0
  
  call VecWAXPY(field%tran_dxx,-1.d0,field%tran_xx,field%tran_yy,ierr)
  
  call VecStrideNorm(field%tran_dxx,ZERO_INTEGER,NORM_INFINITY,option%dcmax,ierr)
      
end subroutine RTMaxChange

! ************************************************************************** !
!
! RTSetPlotVariables: Adds variables to be printed to list
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine RTSetPlotVariables(realization)
  
  use Realization_module
  use Option_module
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXWORDLENGTH) :: name,  units
  type(output_variable_list_type), pointer :: list
  
  character(len=MAXHEADERLENGTH) :: header
  character(len=MAXSTRINGLENGTH) string
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  PetscInt :: i
  
  option => realization%option
  reaction => realization%reaction
  list => realization%output_option%output_variable_list
  
  if (reaction%print_free_conc_type == PRIMARY_MOLALITY) then
    free_mol_char = 'm'
  else
    free_mol_char = 'M'
  endif
  
  if (reaction%print_tot_conc_type == TOTAL_MOLALITY) then
    tot_mol_char = 'm'
  else
    tot_mol_char = 'M'
  endif
  
  if (reaction%print_secondary_conc_type == SECONDARY_MOLALITY) then
    sec_mol_char = 'm'
  else
    sec_mol_char = 'M'
  endif
  
  if (reaction%print_pH .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%h_ion_id > 0) then
      name = 'pH'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,PH, &
                                   reaction%species_idx%h_ion_id)
    endif
  endif  
  
  if (reaction%print_total_component) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Total ' // trim(reaction%primary_species_names(i))
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     reaction%print_tot_conc_type,i)
      endif
    enddo
  endif  
  
  if (reaction%print_free_ion) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Free ' // trim(reaction%primary_species_names(i)) 
        units = trim(free_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                      reaction%print_free_conc_type,i)
      endif
    enddo
  endif  
  
  if (reaction%print_total_bulk) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Total Bulk ' // trim(reaction%primary_species_names(i))
        units = 'mol/m^3 bulk'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_BULK,i)   
      endif
    enddo
  endif
  
  if (reaction%print_act_coefs) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Gamma ' // trim(reaction%primary_species_names(i))
        units = ''
        call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                     PRIMARY_ACTIVITY_COEF,i) 
      endif
    enddo
  endif
  
  do i=1,reaction%neqcplx
    if (reaction%secondary_species_print(i)) then
      name = trim(reaction%secondary_species_names(i))
      units = trim(sec_mol_char)
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   reaction%print_secondary_conc_type,i) 
    endif
  enddo   

  do i=1,reaction%mineral%nkinmnrl
    if (reaction%mineral%kinmnrl_print(i)) then
      name = trim(reaction%mineral%kinmnrl_names(i)) // ' VF'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_VOLUME_FRACTION,units, &
                                   MINERAL_VOLUME_FRACTION,i)     
    endif
  enddo
  
  do i=1,reaction%mineral%nkinmnrl
    if (reaction%mineral%kinmnrl_print(i)) then
      name = trim(reaction%mineral%kinmnrl_names(i)) // ' Rate'
      units = 'mol/sec'
      call OutputVariableAddToList(list,name,OUTPUT_RATE,units, &
                                   MINERAL_RATE,i)      
    endif
  enddo  
  
  do i=1,reaction%mineral%nmnrl
    if (reaction%mineral%mnrl_print(i)) then
      name = trim(reaction%mineral%kinmnrl_names(i)) // ' SI'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                   MINERAL_SATURATION_INDEX,i)    
    endif
  enddo
  
  do i=1,realization%reaction%surface_complexation%nsrfcplxrxn
    if (reaction%surface_complexation%srfcplxrxn_site_density_print(i)) then
      name = trim(reaction%surface_complexation%srfcplxrxn_site_names(i)) // &
             ' Site Density'
      units = 'mol/m^3 bulk'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   SURFACE_SITE_DENSITY,i)
    endif
  enddo  

  do i=1,realization%reaction%surface_complexation%nsrfcplxrxn
    if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
      name = 'Free ' // &
             trim(reaction%surface_complexation%srfcplxrxn_site_names(i))
      units = 'mol/m^3 bulk'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   SURFACE_CMPLX_FREE,i)
    endif
  enddo
  
  
  do i=1,realization%reaction%surface_complexation%nsrfcplx
    if (reaction%surface_complexation%srfcplx_print(i)) then
      name = reaction%surface_complexation%srfcplx_names(i)
      units = 'mol/m^3 bulk'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   SURFACE_CMPLX,i)
    endif
  enddo

  do i=1,realization%reaction%surface_complexation%nkinsrfcplxrxn
    if (reaction%surface_complexation%srfcplxrxn_site_print(i)) then
      option%io_buffer = 'Printing of kinetic surface complexes needs to be fixed'
      call printErrMsg(option)
      name = 'Free ' // &
             trim(reaction%surface_complexation%srfcplxrxn_site_names(i))
      units = 'mol/m^3 bulk'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   KIN_SURFACE_CMPLX_FREE,i)
    endif
  enddo  
  
  do i=1,realization%reaction%surface_complexation%nkinsrfcplx
    if (reaction%surface_complexation%srfcplx_print(i)) then
      option%io_buffer = 'Printing of kinetic surface complexes needs to be fixed'
      call printErrMsg(option)
      name = reaction%surface_complexation%srfcplx_names(i)
      units = 'mol/m^3 bulk'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   KIN_SURFACE_CMPLX,i)
    endif
  enddo  

  if (associated(reaction%kd_print)) then
    do i=1,reaction%naqcomp
      if (reaction%kd_print(i)) then
      name = trim(reaction%primary_species_names(i)) // ' KD'
      units = '-'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   PRIMARY_KD,i)
      endif
    enddo
  endif
  
  if (associated(reaction%total_sorb_print)) then
    do i=1,reaction%naqcomp
      if (reaction%total_sorb_print(i)) then
        name = 'Total Sorbed ' // trim(reaction%primary_species_names(i))
        units = 'mol/m^3'
       call  OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_SORBED,i)        
      endif
    enddo
  endif
  
  if (associated(reaction%total_sorb_mobile_print)) then
    do i=1,reaction%ncollcomp
      if (reaction%total_sorb_mobile_print(i)) then
        name = 'Total Sorbed Mobile ' // &
               trim(reaction%colloid_species_names(i))
        units = '??'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_SORBED_MOBILE,i)  
      endif
    enddo
  endif  
  
  if (reaction%print_colloid) then
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        name = 'Mobile Colloidal ' // trim(reaction%colloid_names(i)) 
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     COLLOID_MOBILE,i)         
      endif
    enddo
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        name = 'Mobile Colloidal ' // trim(reaction%colloid_names(i))
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     COLLOID_IMMOBILE,i)         
      endif
    enddo
  endif
  
  
  if (reaction%print_age) then
    if (reaction%species_idx%tracer_age_id > 0) then
      name = 'Tracer Age'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                   AGE,reaction%species_idx%tracer_age_id, &
                                   reaction%species_idx%tracer_aq_id)       
    endif
  endif  
  
end subroutine RTSetPlotVariables

! ************************************************************************** !
!
! RTJumpStartKineticSorption: Calculates the concentrations of species sorbing
!                             through kinetic sorption processes based
!                             on equilibrium with the aqueous phase.
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RTJumpStartKineticSorption(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTJumpStartKineticSorptionPatch(realization)
                              
end subroutine RTJumpStartKineticSorption

! ************************************************************************** !
!
! RTJumpStartKineticSorptionPatch: Calculates the concentrations of species 
!                                  sorbing through kinetic sorption processes 
!                                  within a patch based on equilibrium with
!                                  the aqueous phase.
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RTJumpStartKineticSorptionPatch(realization)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction

  PetscInt :: ghosted_id
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  reaction => realization%reaction
  
  ! This subroutine assumes that the auxiliary variables are current!

  if (reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RJumpStartKineticSorption(patch%aux%RT%aux_vars(ghosted_id), &
                                     patch%aux%Global%aux_vars(ghosted_id), &
                                     reaction,option)
    enddo
  endif

end subroutine RTJumpStartKineticSorptionPatch

! ************************************************************************** !
!
! RTCheckpointKineticSorption: Checkpoints expliclity stored sorbed 
!                              concentrations
! author: Glenn Hammond
! date: 08/06/09
!
! ************************************************************************** !
subroutine RTCheckpointKineticSorption(realization,viewer,checkpoint)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscBool :: checkpoint
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, pointer :: vec_p(:)

  PetscBool :: checkpoint_flag(realization%reaction%naqcomp)
  PetscInt :: i, j, irxn, icomp, icplx, ncomp, ncplx, irate, ikinmrrxn
  PetscInt :: local_id
  PetscErrorCode :: ierr
  
  option => realization%option
  reaction => realization%reaction
  field => realization%field
  patch => realization%patch
  
  checkpoint_flag = PETSC_FALSE

  ! Loop over sorption reactions to find the necessary components
  
  do ikinmrrxn = 1, reaction%surface_complexation%nkinmrsrfcplxrxn
    irxn = reaction%surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn(ikinmrrxn)
    ncplx = reaction%surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do j = 1, ncplx
      icplx = reaction%surface_complexation%srfcplxrxn_to_complex(j,irxn)
      ncomp = reaction%surface_complexation%srfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%surface_complexation%srfcplxspecid(i,icplx)
        checkpoint_flag(icomp) = PETSC_TRUE
      enddo
    enddo
  enddo

  rt_auxvars => patch%aux%RT%aux_vars
  grid => patch%grid
  do icomp = 1, reaction%naqcomp
    if (checkpoint_flag(icomp)) then
      do irxn = 1, reaction%surface_complexation%nkinmrsrfcplxrxn
        do irate = 1, reaction%surface_complexation%kinmr_nrate(irxn)
          if (checkpoint) then
            call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
            do local_id = 1, grid%nlmax
              vec_p(local_id) = &
                rt_auxvars(grid%nL2G(local_id))% &
                  kinmr_total_sorb(icomp,irate,irxn)
            enddo
            call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
            call VecView(field%work,viewer,ierr)
          else
            call VecLoad(field%work,viewer,ierr)
            if (.not.option%no_restart_kinetic_sorption) then
              call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
              do local_id = 1, grid%nlmax
                rt_auxvars(grid%nL2G(local_id))% &
                  kinmr_total_sorb(icomp,irate,irxn) = &
                    vec_p(local_id)
              enddo
              call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
            endif
          endif
        enddo
      enddo
    endif
  enddo

end subroutine RTCheckpointKineticSorption

! ************************************************************************** !
!
! RTExplicitAdvection: Updates advective transport explicitly
! author: Glenn Hammond
! date: 02/03/12
!
! ************************************************************************** !
subroutine RTExplicitAdvection(realization)

  use Realization_module

  type(realization_type) :: realization
  
  call RTExplicitAdvectionPatch(realization)

end subroutine RTExplicitAdvection

! ************************************************************************** !
!
! RTExplicitAdvectionPatch: Updates advective transport explicitly
! author: Glenn Hammond
! date: 02/03/12
!
! ************************************************************************** !
subroutine RTExplicitAdvectionPatch(realization)

  use Realization_module
  use Discretization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(discretization_type), pointer :: discretization
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscInt :: iphase
  PetscInt :: id_up2, id_dn2
  PetscInt :: local_start, local_end, istart, iend
  PetscInt :: ntvddof
  PetscReal :: qsrc, coef_in, coef_out
  PetscReal :: velocity, area, psv_t  
  PetscReal :: flux(realization%reaction%ncomp)
  
  PetscReal :: sum_flux(realization%reaction%ncomp,realization%patch%grid%ngmax)
  
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: tvd_ghosts_p(:)
  PetscReal, pointer :: rhs_coef_p(:)
  PetscReal, pointer :: total_up2(:,:), total_dn2(:,:)
  PetscErrorCode :: ierr
  PetscViewer :: viewer

#ifdef FORTRAN_2003_COMPLIANT  
  procedure (TFluxLimiterDummy), pointer :: TFluxLimitPtr
  
  select case(realization%option%tvd_flux_limiter)
    case(TVD_LIMITER_UPWIND)
      TFluxLimitPtr => TFluxLimitUpwind
    case(TVD_LIMITER_MC)
      TFluxLimitPtr => TFluxLimitMC
    case(TVD_LIMITER_MINMOD)
      TFluxLimitPtr => TFluxLimitMinmod
    case(TVD_LIMITER_SUPERBEE)
      TFluxLimitPtr => TFluxLimitSuperBee
    case(TVD_LIMITER_VAN_LEER)
      TFluxLimitPtr => TFluxLimitVanLeer
    case default
      TFluxLimitPtr => TFluxLimiter
  end select
#endif

  option => realization%option
  field => realization%field
  patch => realization%patch
  discretization => realization%discretization
  reaction => realization%reaction
  grid => patch%grid
!  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  
  ntvddof = patch%aux%RT%rt_parameter%naqcomp
  
  if (realization%option%tvd_flux_limiter /= TVD_LIMITER_UPWIND) then
    allocate(total_up2(option%nphase,ntvddof))
    allocate(total_dn2(option%nphase,ntvddof))
  else
    ! these must be nullifed so that the explicit scheme ignores them
    nullify(total_up2)
    nullify(total_up2)
  endif

  ! load total component concentrations into tran_xx_p.  it will be used
  ! as local storage here and eventually be overwritten upon leaving 
  ! this routine
  call GridVecGetArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  tran_xx_p = 0.d0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    local_end = local_id*ntvddof
    local_start = local_end-ntvddof+1
    do iphase = 1, option%nphase
      tran_xx_p(local_start:local_end) = &
        rt_aux_vars(ghosted_id)%total(:,iphase)
    enddo
  enddo
  call GridVecRestoreArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call VecScatterBegin(discretization%tvd_ghost_scatter,field%tran_xx, &
                       field%tvd_ghosts,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(discretization%tvd_ghost_scatter,field%tran_xx, &
                     field%tvd_ghosts,INSERT_VALUES,SCATTER_FORWARD,ierr)

! Update Boundary Concentrations------------------------------
  call GridVecGetArrayF90(grid,field%tvd_ghosts,tvd_ghosts_p,ierr)
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (associated(cur_connection_set%id_dn2)) then
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 < 0) then
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          tvd_ghosts_p(istart:iend) = rt_aux_vars_bc(sum_connection)%total(1,:)
        endif
      enddo
    endif
    boundary_condition => boundary_condition%next
  enddo  
  call GridVecRestoreArrayF90(grid,field%tvd_ghosts,tvd_ghosts_p,ierr)
#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'tvd_ghosts.out', &
                            viewer,ierr)
  call VecView(field%tvd_ghosts,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  sum_flux = 0.d0
  
  if (reaction%ncoll > 0) then
    option%io_buffer = &
      'Need to add colloidal source/sinks to RTExplicitAdvectionPatch()'
    call printErrMsg(option)
  endif
  if (option%nphase > 1) then
    option%io_buffer = &
      'Need to add multiphase source/sinks to RTExplicitAdvectionPatch()'
    call printErrMsg(option)
  endif
  if (reaction%ncomp /= reaction%naqcomp) then
    option%io_buffer = &
      'Need to account for non-aqueous species to RTExplicitAdvectionPatch()'
    call printErrMsg(option)
  endif
  if (option%compute_mass_balance_new) then  
    option%io_buffer = &
      'Mass balance not yet supported in RTExplicitAdvectionPatch()'
    call printErrMsg(option)
  endif
  
! Interior Flux Terms -----------------------------------
  call GridVecGetArrayF90(grid,field%tvd_ghosts,tvd_ghosts_p,ierr)
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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
        
      if (associated(cur_connection_set%id_dn2)) then
        id_up2 = cur_connection_set%id_up2(iconn)
        if (id_up2 > 0) then
          total_up2 = rt_aux_vars(id_up2)%total
        else
          iend = abs(id_up2)*ntvddof
          istart = iend-ntvddof+1
          total_up2(1,:) = tvd_ghosts_p(istart:iend)
        endif
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 > 0) then
          total_dn2 = rt_aux_vars(id_dn2)%total
        else
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          total_dn2(1,:) = tvd_ghosts_p(istart:iend)
        endif
      endif
      call TFluxTVD(patch%aux%RT%rt_parameter, &
                    patch%internal_velocities(:,sum_connection), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    total_up2, &
                    rt_aux_vars(ghosted_id_up), &
                    rt_aux_vars(ghosted_id_dn), &
                    total_dn2, &
#ifdef FORTRAN_2003_COMPLIANT  
                    TFluxLimitPtr, &
#endif      
                    option,flux)
          
      ! contribution upwind
      sum_flux(:,ghosted_id_up) = sum_flux(:,ghosted_id_up) - flux
        
      ! contribution downwind
      sum_flux(:,ghosted_id_dn) = sum_flux(:,ghosted_id_dn) + flux
          
    enddo ! iconn
    cur_connection_set => cur_connection_set%next
  enddo
  call GridVecRestoreArrayF90(grid,field%tvd_ghosts,tvd_ghosts_p,ierr)
    
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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(cur_connection_set%id_dn2)) then
        total_up2 = rt_aux_vars_bc(sum_connection)%total
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 > 0) then
          total_dn2 = rt_aux_vars(id_dn2)%total
        else
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          total_dn2(1,:) = tvd_ghosts_p(istart:iend)
        endif
      endif
      call TFluxTVD(patch%aux%RT%rt_parameter, &
                    patch%boundary_velocities(:,sum_connection), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    total_up2, &
                    rt_aux_vars_bc(sum_connection), &
                    rt_aux_vars(ghosted_id), &
                    total_dn2, &
#ifdef FORTRAN_2003_COMPLIANT  
                    TFluxLimitPtr, &
#endif      
                    option,flux)

      ! contribution downwind
      sum_flux(:,ghosted_id) = sum_flux(:,ghosted_id) + flux
#if 0      
      
      do iphase = 1, option%nphase
        velocity = patch%boundary_velocities(iphase,sum_connection)
        area = cur_connection_set%area(iconn)
        
        if (velocity > 0.d0) then  ! inflow
          flux = velocity*area* &
                  rt_aux_vars_bc(sum_connection)%total(:,iphase)
        else  ! outflow
          flux = velocity*area* &
                  rt_aux_vars(ghosted_id)%total(:,iphase)
        endif
          
        ! contribution downwind
        sum_flux(:,ghosted_id) = sum_flux(:,ghosted_id) + flux
          
      enddo ! iphase
#endif        
     
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
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

      do iphase = 1, option%nphase
        qsrc = patch%ss_fluid_fluxes(iphase,sum_connection)
        call TSrcSinkCoef(option,qsrc,source_sink%tran_condition%itype, &
                          coef_in,coef_out)
        flux = coef_in*rt_aux_vars(ghosted_id)%total(:,iphase) + &
               coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                          rt_auxvar%total(:,iphase)
        !geh: TSrcSinkCoef() unit are in L/s.
         sum_flux(:,ghosted_id) = sum_flux(:,ghosted_id) + flux
      enddo
    enddo
    source_sink => source_sink%next
  enddo
  
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)

  
  ! update concentration
  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    local_end = local_id*ntvddof
    local_start = local_end-ntvddof+1
!    do iphase = 1, option%nphase
      ! psv_t must have same units [mol/sec] and be consistent with rhs_coef_p
      ! in RTUpdateRHSCoefsPatch()
      psv_t = porosity_loc_p(ghosted_id)* &
              global_aux_vars(ghosted_id)%sat(iphase)* &
              1000.d0* &
              volume_p(local_id)/option%tran_dt
      !geh: clearly dangerous that I reload into total, but I am going to do it!
      tran_xx_p(local_start:local_end) = &
        ((rhs_coef_p(local_id)*rt_aux_vars(ghosted_id)%total(:,iphase)) + &
         sum_flux(:,ghosted_id)) / psv_t
!    enddo
  enddo
  
  if (associated(total_up2)) then
    deallocate(total_up2)
    nullify(total_up2)
    deallocate(total_dn2)
    nullify(total_dn2)
  endif
  
  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)
  
end subroutine RTExplicitAdvectionPatch

! ************************************************************************** !
!
! RTAppendToHeader: Appends formatted strings to header string
! author: Glenn Hammond
! date: 10/27/11
!
! ************************************************************************** !
subroutine RTAppendToHeader(header,variable_string,cell_string,icolumn)

  character(len=MAXHEADERLENGTH) :: header
  character(len=*) :: variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXSTRINGLENGTH) :: variable_string_adj
  character(len=MAXWORDLENGTH) :: column_string
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: len_cell_string

  variable_string_adj = variable_string
  !geh: Shift to left.  Cannot perform on same string since len=*
  variable_string_adj = adjustl(variable_string_adj)
  
  if (icolumn > 0) then
    icolumn = icolumn + 1
    write(column_string,'(i4,''-'')') icolumn
    column_string = trim(adjustl(column_string))
  else
    column_string = ''
  endif

  !geh: this is all to remove the lousy spaces
  len_cell_string = len_trim(cell_string) 

  if (len_cell_string > 0) then
    write(string,'('',"'',a,a,'' '',a,''"'')') trim(column_string), &
          trim(variable_string_adj), trim(cell_string)
  else
    write(string,'('',"'',a,a,''"'')') trim(column_string), &
          trim(variable_string_adj)
  endif
  header = trim(header) // trim(string)

end subroutine RTAppendToHeader

! ************************************************************************** !
!
! RTSecondaryTransportMulti:  Calculates the source term contribution due to 
! secondary continuum in the primary continuum residual for multicomponent
! system assuming only aqueous reaction
! author: Satish Karra
! date: 1/31/13
!
! ************************************************************************** !
subroutine RTSecondaryTransportMulti(sec_transport_vars,aux_var, &
                                     global_aux_var, &
                                     reaction,diffusion_coefficient, &
                                     porosity,option,res_transport)
                               
                            
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module
  use blksolv_module
  use Utility_module
  

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: aux_var
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type) :: global_aux_var
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells-1)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells-1)
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: identity(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: b_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: sec_jac(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: inv_D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: conc_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: conc_current_M(reaction%naqcomp)
  PetscReal :: total_current_M(reaction%naqcomp)
  PetscReal :: res_transport(reaction%naqcomp)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)

  
  PetscInt :: i, j, k, n
  PetscInt :: ngcells, ncomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal, parameter :: rgas = 8.3144621d-3
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  conc_upd = sec_transport_vars%updated_conc
  ncomp = reaction%naqcomp

  do j = 1, ncomp
    do i = 1, ngcells
      conc_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j)* &
                       global_aux_var%den_kg(1)*1.d-3 
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
    enddo
  enddo
  
  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  res = 0.d0
  rhs = 0.d0
  D_M = 0.d0
  identity = 0.d0
  b_M = 0.d0
  inv_D_M = 0.d0
  total_upd = 0.d0 ! Need to calculate this from conc_upd
  total_current_M = 0.d0
  
  
  total_primary_node = aux_var%total(:,1)                             ! in mol/L 
  pordt = porosity/option%tran_dt
  pordiff = porosity*diffusion_coefficient            
              
!================ Calculate the secondary residual =============================        

  do j = 1, ncomp
      
    ! Accumulation
    do i = 1, ngcells
      n = j + (i-1)*ncomp
      res(n) = pordt*(total_upd(j,i) - total_prev(j,i))*vol(i)    ! in mol/L*m3/s
    enddo
  
    ! Flux terms
    do i = 2, ngcells - 1
      n = j + (i-1)*ncomp
      res(n) = res(n) - pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (total_upd(j,i+1) - total_upd(j,i))
      res(n) = res(n) + pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (total_upd(j,i) - total_upd(j,i-1))                      
    enddo
         
              
    ! Apply boundary conditions
    ! Inner boundary
    res(j) = res(j) - pordiff*area(1)/(dm_minus(2) + dm_plus(1))* &
                      (total_upd(j,2) - total_upd(j,1))
                                      
    ! Outer boundary
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) - &
                               pordiff*area(ngcells)/dm_plus(ngcells)* &
                               (total_primary_node(j) - total_upd(j,ngcells))
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) + &
                               pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                               + dm_plus(ngcells-1))*(total_upd(j,ngcells) - &
                               total_upd(j,ngcells-1))  
                               
  enddo
                         
  res = res*1.d3 ! Convert mol/L*m3/s to mol/s                                                    
                                                                                                          
!================ Calculate the secondary jacobian =============================        


  do j = 1, ncomp
    do k = 1, ncomp  
      if (k == j) then
        ! Accumulation
        do i = 1, ngcells 
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + pordt*vol(i)
        enddo
  
        ! Flux terms
        do i = 2, ngcells-1
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + &
                              pordiff*area(i)/(dm_minus(i+1) + dm_plus(i)) + &
                              pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_left(j,k,i-1) = coeff_left(j,k,i-1) - &
                              pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_right(j,k,i) = coeff_right(j,k,i) - &
                               pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))
        enddo
  
  
        ! Apply boundary conditions
        ! Inner boundary
        coeff_diag(j,k,1) = coeff_diag(j,k,1) + &
                            pordiff*area(1)/(dm_minus(2) + dm_plus(1))
                   
        coeff_right(j,k,1) = coeff_right(j,k,1) - &
                             pordiff*area(1)/(dm_minus(2) + dm_plus(1))
  
        ! Outer boundary -- closest to primary node
        coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                  pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                                  + dm_plus(ngcells-1)) + &
                                  pordiff*area(ngcells)/dm_plus(ngcells)
        coeff_left(j,k,ngcells-1) = coeff_left(j,k,ngcells-1) - &
                                  pordiff*area(ngcells-1)/(dm_minus(ngcells) + &
                                  dm_plus(ngcells-1)) 

      endif
    enddo    
  enddo
  
  ! Convert mol/L*m3/s to mol/s
  coeff_right = coeff_right*1.d3
  coeff_left = coeff_left*1.d3
  coeff_diag = coeff_diag*1.d3
         
!===============================================================================        
                        
  rhs = -res                 
 
  call bl3dfac(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot)
  
  call bl3dsolf(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot,1,rhs)
    
  ! Set the values of D_M matrix and create identity matrix of size ncomp x ncomp
  do i = 1, ncomp
    do j = 1, ncomp
      D_M(i,j) = coeff_diag(i,j,ngcells)
      if (j == i) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo

  ! Find the inverse of D_M
  call ludcmp(D_M,ncomp,indx,d) 
  do j = 1, ncomp
    call lubksb(D_M,ncomp,indx,identity(1,j))
  enddo  
  inv_D_M = identity 
  
  ! Update the secondary concentrations
  do i = 1, ncomp
    conc_current_M(i) = conc_upd(i,ngcells) + rhs(i+(ngcells-1)*ncomp)
  enddo
  
  ! Update the secondary continuum totals at the outer matrix node
  rt_auxvar => sec_transport_vars%sec_rt_auxvar(ngcells)
  rt_auxvar%pri_molal(1:ncomp) = conc_current_M/(global_aux_var%den_kg(1)*1.d-3) ! in mol/kg
  call RTotal(rt_auxvar,global_aux_var,reaction,option)
  
  total_current_M = rt_auxvar%total(:,1)
  
  b_m = pordiff/dm_plus(ngcells)*area(ngcells)*inv_D_M
  
  ! Reset identity matrix
  do i = 1, ncomp
    do j = 1, ncomp
      if (j == i) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo
  
  ! Calculate the coupling term
  res_transport = pordiff/dm_plus(ngcells)*area_fm* &
                  (conc_current_M - total_primary_node)
                  
  ! Calculate the jacobian contribution due to coupling term
  sec_jac = area_fm*pordiff/dm_plus(ngcells)*(b_m - identity)
    
  ! Store the contribution to the primary jacobian term
  sec_transport_vars%sec_jac = sec_jac 
  sec_transport_vars%sec_jac_update = PETSC_TRUE
  
  ! Store the coefficients from LU decomposition of the block tridiagonal
  ! sytem. These will be called later to perform backsolve to the get the
  ! updated secondary continuum concentrations at the end of the timestep
  sec_transport_vars%cxm = coeff_left
  sec_transport_vars%cxp = coeff_right
  sec_transport_vars%cdl = coeff_diag
  
  ! Store the solution of the forward solve
  sec_transport_vars%r = rhs

end subroutine RTSecondaryTransportMulti


! ************************************************************************** !
!
! RTDestroy: Deallocates variables associated with Reactive Transport
! author: Glenn Hammond
! date: 02/03/09
!
! ************************************************************************** !
subroutine RTDestroy(realization)

  use Realization_module
  use Patch_module
  use Option_module

  type(realization_type) :: realization
  
#ifdef OS_STATISTICS
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  PetscReal :: temp_real_in(3), temp_real_out(3)
  PetscReal :: call_count
  PetscReal :: sum_newton_iterations
  PetscReal :: ave_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_in_a_cell
  PetscReal :: max_newton_iterations_on_a_core
  PetscReal :: min_newton_iterations_on_a_core
  
  PetscReal :: sum, ave, var, value
  PetscInt :: irank
  PetscReal, allocatable :: tot_newton_iterations(:)
  
  option => realization%option
  call_count = 0.d0
  sum_newton_iterations = 0.d0
  max_newton_iterations_in_a_cell = -99999999
  max_newton_iterations_on_a_core = -99999999.d0
  min_newton_iterations_on_a_core = 99999999.d0
  
#endif  

#ifdef OS_STATISTICS
  call_count = call_count + &
    cur_patch%aux%RT%rt_parameter%sum_newton_call_count
  sum_newton_iterations = sum_newton_iterations + &
    cur_patch%aux%RT%rt_parameter%sum_newton_iterations
  if (cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations > &
      max_newton_iterations_in_a_cell) then
    max_newton_iterations_in_a_cell = &
      cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations
  endif
#endif
  call RTDestroyPatch(realization)

#ifdef OS_STATISTICS
  if (option%reactive_transport_coupling == OPERATOR_SPLIT) then
    temp_real_in(1) = call_count
    temp_real_in(2) = sum_newton_iterations
    call MPI_Allreduce(temp_real_in,temp_real_out,TWO_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    ave_newton_iterations_in_a_cell = temp_real_out(2)/temp_real_out(1)

    temp_real_in(1) = dble(max_newton_iterations_in_a_cell)
    temp_real_in(2) = sum_newton_iterations
    temp_real_in(3) = -sum_newton_iterations
    call MPI_Allreduce(temp_real_in,temp_real_out,THREE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    max_newton_iterations_in_a_cell = int(temp_real_out(1)+1.d-4)
    max_newton_iterations_on_a_core = temp_real_out(2)
    min_newton_iterations_on_a_core = -temp_real_out(3)

    ! Now let's compute the variance!
    call OptionMeanVariance(sum_newton_iterations,ave,var,PETSC_TRUE,option)
  
    if (option%print_screen_flag) then
      write(*, '(/,/" OS Reaction Statistics (Overall): ",/, &
               & "       Ave Newton Its / Cell: ",1pe12.4,/, &
               & "       Max Newton Its / Cell: ",i4,/, &
               & "       Max Newton Its / Core: ",1pe12.4,/, &
               & "       Min Newton Its / Core: ",1pe12.4,/, &
               & "       Ave Newton Its / Core: ",1pe12.4,/, &
               & "   Std Dev Newton Its / Core: ",1pe12.4,/)') &
                 ave_newton_iterations_in_a_cell, &
                 max_newton_iterations_in_a_cell, &
                 max_newton_iterations_on_a_core, &
                 min_newton_iterations_on_a_core, &
                 ave, &
                 sqrt(var)

    endif

    if (option%print_file_flag) then
      write(option%fid_out, '(/,/" OS Reaction Statistics (Overall): ",/, &
               & "       Ave Newton Its / Cell: ",1pe12.4,/, &
               & "       Max Newton Its / Cell: ",i4,/, &
               & "       Max Newton Its / Core: ",1pe12.4,/, &
               & "       Min Newton Its / Core: ",1pe12.4,/, &
               & "       Ave Newton Its / Core: ",1pe12.4,/, &
               & "   Std Dev Newton Its / Core: ",1pe12.4,/)') &
                 ave_newton_iterations_in_a_cell, &
                 max_newton_iterations_in_a_cell, &
                 max_newton_iterations_on_a_core, &
                 min_newton_iterations_on_a_core, &
                 ave, &
                 sqrt(var)
    endif
  endif

#endif 


end subroutine RTDestroy

! ************************************************************************** !
!
! RTDestroyPatch: Deallocates variables associated with Reactive Transport
! author: Glenn Hammond
! date: 02/03/09
!
! ************************************************************************** !
subroutine RTDestroyPatch(realization)

  use Realization_module

  implicit none

  type(realization_type) :: realization
  
  ! taken care of in auxiliary.F90
  
end subroutine RTDestroyPatch

end module Reactive_Transport_module
