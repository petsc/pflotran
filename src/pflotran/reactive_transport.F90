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
  
  public :: RTTimeCut, RTSetup, RTMaxChange, RTUpdateSolution, RTResidual, &
            RTJacobian, RTInitializeTimestep, RTGetTecplotHeader, &
            RTUpdateAuxVars, RTComputeMassBalance, RTDestroy, RTCheckUpdate
  
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
  call VecCopy(field%tran_xx,field%tran_yy,ierr)
  
  ! set densities and saturations to t
  if (realization%option%nflowdof > 0) then
    call GlobalUpdateDenAndSat(realization,realization%option%tran_weight_t0)
  endif
  
  call RTInitializeTimestep(realization)  
  
  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalUpdateDenAndSat(realization,realization%option%tran_weight_t1)
  endif
 
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
      cur_patch%reaction => realization%reaction
      call RTSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  use Fluid_module
  use Material_module
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(fluid_property_type), pointer :: cur_fluid_property

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: iphase
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction

  patch%aux%RT => RTAuxCreate(option)
    
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
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  option%iflag = 1 ! enable allocation of mass_balance array 
  allocate(patch%aux%RT%aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call RTAuxVarInit(patch%aux%RT%aux_vars_bc(iconn),reaction,option)
  enddo
  patch%aux%RT%num_aux_bc = sum_connection
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
    cur_fluid_property => cur_fluid_property%next
  enddo
  
  if (associated(realization%material_properties)) then
    patch%aux%RT%rt_parameter%dispersivity = &
      realization%material_properties%longitudinal_dispersivity
  endif

end subroutine RTSetupPatch

! ************************************************************************** !
!
! RTCheckUpdate: In the case of the log formulation, ensures that the update 
!                vector does not exceed a prescribed tolerance
! author: Glenn Hammond
! date: 03/16/09
!
! ************************************************************************** !
subroutine RTCheckUpdate(snes_,C,dC,realization,changed,ierr)
 
  use Realization_module
  use Level_module
  use Patch_module
 
  implicit none
  
  SNES :: snes_
  Vec :: C
  Vec :: dC
  type(realization_type) :: realization
  PetscTruth :: changed
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
      call RTCheckUpdatePatch(snes_,C,dC,realization,changed,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTCheckUpdate

! ************************************************************************** !
!
! RTCheckUpdatePatch: In the case of the log formulation, ensures that the 
!                     update vector does not exceed a prescribed tolerance
! author: Glenn Hammond
! date: 03/16/09
!
! ************************************************************************** !
subroutine RTCheckUpdatePatch(snes_,C,dC,realization,changed,ierr)

  use Realization_module
  use Grid_module
 
  implicit none
  
  SNES :: snes_
  Vec :: C
  Vec :: dC
  type(realization_type) :: realization
  PetscTruth :: changed
  
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
    call MPI_AllReduce(ratio,min_ratio,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       PETSC_COMM_WORLD,ierr)
                       
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
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%ntrandof, &
                            realization%option%nphase)
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
    
  mass_balance = 0.d0
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTComputeMassBalancePatch(realization,mass_balance)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  PetscInt :: i, icomp, imnrl, ncomp, irate

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
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    do iphase = 1, option%nphase
      mass_balance(:,iphase) = mass_balance(:,iphase) + &
        rt_aux_vars(ghosted_id)%total(:,iphase) * &
        global_aux_vars(ghosted_id)%sat(iphase) * &
        porosity_loc_p(ghosted_id) * &
        volume_p(local_id)*1000.d0
        
      if (iphase == 1) then
      ! add contribution of equilibrium sorption
        if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
          mass_balance(:,iphase) = mass_balance(:,iphase) + &
            rt_aux_vars(ghosted_id)%total_sorb_eq(:) * volume_p(local_id)
        endif


      ! add contribution of kinetic multirate sorption
        if (reaction%kinmr_nrate > 0) then
          do irate = 1, reaction%kinmr_nrate
            mass_balance(:,iphase) = mass_balance(:,iphase) + &
              rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate) * &
              volume_p(local_id)
          enddo
        endif
               
      ! add contribution from mineral volume fractions
        if (reaction%nkinmnrl > 0) then
          do imnrl = 1, reaction%nkinmnrl
            ncomp = reaction%kinmnrlspecid(0,imnrl)
            do i = 1, ncomp
              icomp = reaction%kinmnrlspecid(i,imnrl)
              mass_balance(icomp,iphase) = mass_balance(icomp,iphase) &
              + reaction%kinmnrlstoich(i,imnrl)                  &
              * rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl)      &
              * volume_p(local_id) &
              / reaction%kinmnrl_molar_vol(imnrl)
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

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%RT%num_aux
    patch%aux%RT%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  do iconn = 1, patch%aux%RT%num_aux_bc
    rt_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
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

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc

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
      call RTInitializeTimestepPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  use Field_module
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy,ierr)
  
  ! update mineral volume fractions
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
 
  implicit none

  type(realization_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  PetscInt :: ghosted_id, imnrl
  PetscInt :: irate
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

  call RTUpdateAuxVarsPatch(realization,PETSC_FALSE,PETSC_FALSE)
  
  if (.not.option%init_stage) then
    ! update mineral volume fractions
    if (reaction%nkinmnrl > 0) then
      do ghosted_id = 1, grid%ngmax
        do imnrl = 1, reaction%nkinmnrl
          ! rate = mol/m^3/sec
          ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
          !                                mol_vol (m^3 mnrl/mol mnrl)
          rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = &
            rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) + &
            rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* &
            reaction%kinmnrl_molar_vol(imnrl)* &
            option%tran_dt
          if (rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) < 0.d0) &
            rt_aux_vars(ghosted_id)%mnrl_volfrac(imnrl) = 0.d0
        enddo
      enddo
    endif
    
    ! update multirate sorption concentrations 
  ! WARNING: below assumes site concentration multiplicative factor
    if (reaction%kinmr_nrate > 0) then 
      do ghosted_id = 1, grid%ngmax 
        do irate = 1, reaction%kinmr_nrate 
          kdt = reaction%kinmr_rate(irate) * option%tran_dt 
          one_plus_kdt = 1.d0 + kdt 
          k_over_one_plus_kdt = reaction%kinmr_rate(irate)/one_plus_kdt 
          rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate) = & 
            (rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate) + & 
            kdt * reaction%kinmr_frac(irate) * rt_aux_vars(ghosted_id)%total_sorb_eq)/one_plus_kdt 
        enddo 
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

  implicit none
  
  type(realization_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: xx_p(:), porosity_loc_p(:), tor_loc_p(:), &
                        volume_p(:), accum_p(:), density_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid
  reaction => realization%reaction

  ! cannot use tran_xx_loc vector here as it has not yet been updated.
  call GridVecGetArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)

! Do not use RTUpdateAuxVarsPatch() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*reaction%ncomp
    istart = iend-reaction%ncomp+1

    rt_aux_vars(ghosted_id)%pri_molal = xx_p(istart:iend)
    ! DO NOTE RECOMPUTE THE ACTIVITY COEFFICIENTS BEFORE COMPUTING THE
    ! FIXED PORTION OF THE ACCUMULATION TERM - geh
    call RTAuxVarCompute(rt_aux_vars(ghosted_id), &
                         global_aux_vars(ghosted_id), &
                         reaction,option)
    call RTAccumulation(rt_aux_vars(ghosted_id), &
                        global_aux_vars(ghosted_id), &
                        porosity_loc_p(ghosted_id), &
                        volume_p(local_id), &
                        reaction,option,accum_p(istart:iend)) 
    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
      call RAccumulationSorb(rt_aux_vars(ghosted_id), &
                             global_aux_vars(ghosted_id), &
                             volume_p(local_id), &
                             reaction,option,accum_p(istart:iend))
    endif
  enddo

  call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)

end subroutine RTUpdateFixedAccumulationPatch

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
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
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
! RTAccumulation: Computes aqueous portion of the accumulation term in 
!                 residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulation(rt_aux_var,global_aux_var,por,vol,reaction,option,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
  PetscInt :: iphase
  PetscReal :: psv_t
  PetscReal :: v_t
  
  iphase = 1
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  psv_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
  Res(:) = psv_t*rt_aux_var%total(:,iphase) 

! Add in multiphase, clu 12/29/08
#if 1  
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit

! super critical CO2 phase
    if (iphase == 2 ) then
      psv_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
      Res(:) = Res(:) + psv_t*rt_aux_var%total(:,iphase) 
      ! should sum over gas component only need more implementations
    endif 
! add code for other phases here
  enddo

#endif
end subroutine RTAccumulation

! ************************************************************************** !
!
! RTAccumulationDerivative: Computes derivative of aqueous portion of the 
!                           accumulation term in residual function 
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTAccumulationDerivative(rt_aux_var,global_aux_var, &
                                    por,vol,reaction,option,J)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var  
  PetscReal :: por, vol
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp, iphase
  PetscReal :: psvd_t, v_t
  
  iphase = 1
  ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  !         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  if (associated(rt_aux_var%dtotal)) then ! units of dtotal = kg water/L water
    psvd_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
    J = rt_aux_var%dtotal(:,:,iphase)*psvd_t
  else
    J = 0.d0
    psvd_t = por*global_aux_var%sat(iphase)* &
             global_aux_var%den_kg(iphase)*vol/option%tran_dt ! units of den = kg water/m^3 water
    do icomp=1,reaction%ncomp
      J(icomp,icomp) = psvd_t
    enddo
  endif

! Add in multiphase, clu 12/29/08
#if 1  
 do
   iphase = iphase +1 
   if (iphase > option%nphase) exit
! super critical CO2 phase
   if (iphase ==2 ) then
     if (associated(rt_aux_var%dtotal)) then
       psvd_t = por*global_aux_var%sat(iphase)*1000.d0*vol/option%tran_dt  
       J = J + rt_aux_var%dtotal(:,:,iphase)*psvd_t
     else
       J = 0.d0
       psvd_t = por*global_aux_var%sat(iphase)* &
          global_aux_var%den_kg(iphase)*vol/option%tran_dt ! units of den = kg water/m^3 water
       do icomp=1,reaction%ncomp
         J(icomp,icomp) = J(icomp,icomp) + psvd_t
       enddo
     endif   
   endif
 enddo
#endif     
end subroutine RTAccumulationDerivative

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
  use Level_module
  use Discretization_module
  use Option_module
  use Grid_module

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
  PetscReal, pointer :: xx_p(:), log_xx_p(:)
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscViewer :: viewer  
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  if (realization%reaction%use_log_formulation) then
    ! have to convert the log concentration to non-log form
    cur_level => realization%level_list%first
    do
      if (.not.associated(cur_level)) exit
      cur_patch => cur_level%patch_list%first
      do
        if (.not.associated(cur_patch)) exit
        call GridVecGetArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
        call GridVecGetArrayF90(cur_patch%grid,xx,log_xx_p,ierr)
        xx_p(:) = exp(log_xx_p(:))
        call GridVecRestoreArrayF90(cur_patch%grid,field%tran_xx,xx_p,ierr)
        call GridVecRestoreArrayF90(cur_patch%grid,xx,log_xx_p,ierr)  
        cur_patch => cur_patch%next
      enddo
      cur_level => cur_level%next
    enddo
    call DiscretizationGlobalToLocal(discretization,field%tran_xx,field%tran_xx_loc,NTRANDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc,NTRANDOF)
  endif

  ! pass #1 for internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTResidualPatch1(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! now coarsen all face fluxes in case we are using SAMRAI to 
  ! ensure consistent fluxes at coarse-fine interfaces
  if(option%use_samr) then
     call SAMRCoarsenFaceFluxes(discretization%amrgrid%p_application, field%tran_face_fluxes, ierr)

     cur_level => realization%level_list%first
     do
        if (.not.associated(cur_level)) exit
        cur_patch => cur_level%patch_list%first
        do
           if (.not.associated(cur_patch)) exit
           realization%patch => cur_patch
           call RTResidualFluxContribPatch(r,realization,ierr)
           cur_patch => cur_patch%next
        enddo
        cur_level => cur_level%next
     enddo
  endif

  ! pass #2 for everything else
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTResidualPatch2(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif

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
  
end subroutine RTResidual

! ************************************************************************** !
!
! RichardsResidualfuxContribsPatch: should be called only for SAMR
! author: Bobby Philip
! date: 02/17/09
!
! ************************************************************************** !
subroutine RTResidualFluxContribPatch(r,realization,ierr)
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Debug_module
  
  implicit none

  Vec, intent(out) :: r
  type(realization_type) :: realization

  PetscErrorCode :: ierr

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: face_fluxes_p(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscInt :: axis, nlx, nly, nlz
  PetscInt :: iconn, i, j, k
  PetscInt :: xup_id, xdn_id, yup_id, ydn_id, zup_id, zdn_id
  PetscInt :: su1, eu1, sd1, ed1
  PetscInt :: su2, eu2, sd2, ed2
  PetscInt :: su3, eu3, sd3, ed3
  PetscInt :: istart, iend

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction

! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)

  do axis=0,2  
     call GridVecGetArrayF90(grid,axis,field%tran_face_fluxes, fluxes(axis)%flux_p, ierr)  
  enddo

  nlx = grid%structured_grid%nlx  
  nly = grid%structured_grid%nly  
  nlz = grid%structured_grid%nlz 
  
  iconn=1
  do k=1,nlz
     do j=1,nly
        do i=1,nlx
           xup_id = ((k-1)*nly+j-1)*(nlx+1)+i
           xdn_id = xup_id+1
           yup_id = ((k-1)*(nly+1)+(j-1))*nlx+i
           ydn_id = yup_id+nlx
           zup_id = ((k-1)*nly+(j-1))*nlx+i
           zdn_id = zup_id+nlx*nly

           su1    = xup_id*reaction%ncomp
           eu1    = su1+reaction%ncomp-1
           sd1    = xdn_id*reaction%ncomp
           ed1    = sd1+reaction%ncomp-1

           su2    = yup_id*reaction%ncomp
           eu2    = su1+reaction%ncomp-1
           sd2    = ydn_id*reaction%ncomp
           ed2    = sd1+reaction%ncomp-1

           su3    = zup_id*reaction%ncomp
           eu3    = su1+reaction%ncomp-1
           sd3    = zdn_id*reaction%ncomp
           ed3    = sd1+reaction%ncomp-1

           istart=iconn
           iend  = iconn+reaction%ncomp-1
           r_p(istart:iend) = r_p(istart:iend) &
                             +fluxes(0)%flux_p(sd1:ed1)-fluxes(0)%flux_p(su1:eu1) &
                             +fluxes(1)%flux_p(sd2:ed2)-fluxes(1)%flux_p(su2:eu2) &
                             +fluxes(2)%flux_p(sd3:ed3)-fluxes(2)%flux_p(su3:eu3) 

           iconn=iend+1
        enddo
     enddo
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
!!$
!!$  do axis=0,2  
!!$     call GridVecRestoreArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
!!$  enddo

end subroutine RTResidualFluxContribPatch


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
  
  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
     implicit none
     
#include "finclude/petsc.h"
     
     PetscFortranAddr :: p_patch
     PetscInt :: axis,dim
   end function samr_patch_at_bc
  end interface

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
  PetscReal :: Res(realization%reaction%ncomp)
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn

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
  
  if (.not.patch%aux%RT%aux_vars_up_to_date .and. &
      reaction%act_coef_update_frequency == ACT_COEF_FREQUENCY_NEWTON_ITER) then
    call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_TRUE)
  else
    call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_FALSE)
  endif
  patch%aux%RT%aux_vars_up_to_date = PETSC_FALSE 
  
  if (option%compute_mass_balance_new) then
    call RTZeroMassBalanceDeltaPatch(realization)
  endif
  
  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,r, r_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)

  r_p = 0.d0

  if (option%use_samr) then
     do axis=0,2  
        call GridVecGetArrayF90(grid,axis,field%tran_face_fluxes, fluxes(axis)%flux_p, ierr)  
     enddo

    nlx = grid%structured_grid%nlx  
    nly = grid%structured_grid%nly  
    nlz = grid%structured_grid%nlz 

    ngx = grid%structured_grid%ngx   
    ngxy = grid%structured_grid%ngxy

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 0, 0)==1) nlx = nlx-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 0, 1)==1) nlx = nlx-1
    
    max_x_conn = (nlx+1)*nly*nlz
    ! reinitialize nlx
    nlx = grid%structured_grid%nlx  

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 1, 0)==1) nly = nly-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, 1, 1)==1) nly = nly-1
    
    max_y_conn = max_x_conn + nlx*(nly+1)*nlz

    ! reinitialize nly
    nly = grid%structured_grid%nly  

    ! Interior Advective Flux Terms -----------------------------------
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

          call TFluxAdv(rt_aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
               rt_aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
               cur_connection_set%area(iconn),rt_parameter,option, &
               patch%internal_velocities(:,iconn),Res)
          
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
          
       enddo
       cur_connection_set => cur_connection_set%next
    enddo

    ! Interior Diffusive Flux Terms -----------------------------------
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
          dist_up = distance*fraction_upwind
          dist_dn = distance-dist_up ! should avoid truncation error
          
          call TFluxDiff(rt_aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
               porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
               dist_up, &
               rt_aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
               porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
               dist_up, &
               cur_connection_set%area(iconn),rt_parameter,option, &
               patch%internal_velocities(:,iconn),Res)
          
          if (sum_connection <= max_x_conn) then
             direction = 0
             if(mod(mod(ghosted_id_dn,ngxy),ngx).eq.0) then
                flux_id = ((ghosted_id_dn/ngxy)-1)*(nlx+1)*nly + &
                     ((mod(ghosted_id_dn,ngxy))/ngx-1)*(nlx+1)
             else
                flux_id = ((ghosted_id_dn/ngxy)-1)*(nlx+1)*nly + &
                     ((mod(ghosted_id_dn,ngxy))/ngx-1)*(nlx+1)+ &
                     mod(mod(ghosted_id_dn,ngxy),ngx)-1
             endif
             
          else if (sum_connection <= max_y_conn) then
             direction = 1
             flux_id = ((ghosted_id_dn/ngxy)-1)*nlx*(nly+1) + &
                  ((mod(ghosted_id_dn,ngxy))/ngx-1)*nlx + &
                  mod(mod(ghosted_id_dn,ngxy),ngx)-1
          else
             direction = 2
             flux_id = ((ghosted_id_dn/ngxy)-1)*nlx*nly &
                  +((mod(ghosted_id_dn,ngxy))/ngx-1)*nlx &
                  +mod(mod(ghosted_id_dn,ngxy),ngx)-1
          endif
          
          iend = flux_id*reaction%ncomp
          istart = iend-reaction%ncomp+1         
          fluxes(direction)%flux_p(istart:iend) = Res(1:reaction%ncomp)
          
          if (option%store_solute_fluxes) then
             patch%internal_fluxes(iphase,1:reaction%ncomp,iconn) = &
                  Res(1:reaction%ncomp)
          endif
       enddo
       cur_connection_set => cur_connection_set%next
    enddo    
  
    ! Advective Boundary Flux Terms -----------------------------------
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
          
          call TBCFluxAdv(boundary_condition%tran_condition%itype, &
               rt_aux_vars_bc(sum_connection), &
               global_aux_vars_bc(sum_connection), &
               rt_aux_vars(ghosted_id), &
               global_aux_vars(ghosted_id), &
               cur_connection_set%area(iconn), &
               rt_parameter,option, &
               patch%boundary_velocities(:,sum_connection),Res)
          
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
          
       enddo
       boundary_condition => boundary_condition%next
    enddo
    
    ! Diffusive Boundary Flux Terms -----------------------------------
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
          
          call TBCFluxDiff(boundary_condition%tran_condition%itype, &
               rt_aux_vars_bc(sum_connection), &
               global_aux_vars_bc(sum_connection), &
               rt_aux_vars(ghosted_id), &
               global_aux_vars(ghosted_id), &
               porosity_loc_p(ghosted_id), &
               tor_loc_p(ghosted_id), &
               cur_connection_set%dist(0,iconn), &
               cur_connection_set%area(iconn), &
               rt_parameter,option, &
               patch%boundary_velocities(:,sum_connection),Res)
          
          direction =  (boundary_condition%region%faces(iconn)-1)/2
          
          ! the ghosted_id gives the id of the cell. Since the
          ! flux_id is based on the ghosted_id of the downwind
          ! cell this has to be adjusted in the case of the east, 
          ! north and top faces before the flux_id is computed
          select case(boundary_condition%region%faces(iconn)) 
          case(WEST_FACE)
             flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                  ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)+ &
                  mod(mod(ghosted_id,ngxy),ngx)-1
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = Res(1:reaction%ncomp)
          case(EAST_FACE)
             ghosted_id = ghosted_id+1
             flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                  ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = -Res(1:reaction%ncomp)
          case(SOUTH_FACE)
             flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                  ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                  mod(mod(ghosted_id,ngxy),ngx)-1
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = Res(1:reaction%ncomp)
          case(NORTH_FACE)
             ghosted_id = ghosted_id+ngx
             flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                  ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                  mod(mod(ghosted_id,ngxy),ngx)-1
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = -Res(1:reaction%ncomp)
          case(BOTTOM_FACE)
             flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                  +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                  +mod(mod(ghosted_id,ngxy),ngx)-1
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = Res(1:reaction%ncomp)
          case(TOP_FACE)
             ghosted_id = ghosted_id+ngxy
             flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                  +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                  +mod(mod(ghosted_id,ngxy),ngx)-1
             istart = flux_id*reaction%ncomp
             iend = istart+reaction%ncomp-1
             fluxes(direction)%flux_p(istart:iend) = -Res(1:reaction%ncomp)
         end select

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
         
       enddo
       boundary_condition => boundary_condition%next
   enddo

  else
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
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TFlux(rt_aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
                 porosity_loc_p(ghosted_id_up),tor_loc_p(ghosted_id_up), &
                 dist_up, &
                 rt_aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
                 porosity_loc_p(ghosted_id_dn),tor_loc_p(ghosted_id_dn), &
                 dist_up, &
                 cur_connection_set%area(iconn),rt_parameter,option, &
                 patch%internal_velocities(:,iconn),Res)

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

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFlux(boundary_condition%tran_condition%itype, &
                   rt_aux_vars_bc(sum_connection), &
                   global_aux_vars_bc(sum_connection), &
                   rt_aux_vars(ghosted_id), &
                   global_aux_vars(ghosted_id), &
                   porosity_loc_p(ghosted_id), &
                   tor_loc_p(ghosted_id), &
                   cur_connection_set%dist(0,iconn), &
                   cur_connection_set%area(iconn), &
                   rt_parameter,option, &
                   patch%boundary_velocities(:,sum_connection),Res)
 
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
      
    enddo
    boundary_condition => boundary_condition%next
  enddo
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)

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
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscReal :: qsrc, molality
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscTruth :: volumetric
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
  
  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  if (.not.option%steady_state) then
    r_p = r_p - accum_p
#if 1
    ! Accumulation terms ------------------------------------
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      call RTAccumulation(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                          porosity_loc_p(ghosted_id), &
                          volume_p(local_id),reaction,option,Res) 
      r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)
    enddo
  endif
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)
      if (source_sink%flow_condition%rate%itype == &
          VOLUMETRIC_RATE_SS) then
        volumetric = PETSC_TRUE
      else
        volumetric = PETSC_FALSE
      endif
    endif
      
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          ! units should be mol/sec
          Res = -1.d-6* &
                porosity_loc_p(ghosted_id)* &
                global_aux_vars(ghosted_id)%sat(option%liquid_phase)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
                (source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%total(:,iphase) - rt_aux_vars(ghosted_id)%total(:,iphase))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
        case(MASS_RATE_SS)
          Res = -source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%total(:,iphase) ! actually moles/sec
        case default
          if (qsrc > 0) then ! injection
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              Res = -qsrc* &
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%total(:,iphase)*1000.d0
            else
               Res = -qsrc* &
                     source_sink%tran_condition%cur_constraint_coupler% &
                     rt_auxvar%total(:,iphase)/ &
                     global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                     1000.d0
            endif
          else ! extraction
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              Res = -qsrc*rt_aux_vars(ghosted_id)%total(:,iphase)*1000.d0
            else
              Res = -qsrc* &
                    rt_aux_vars(ghosted_id)%total(:,iphase)/ &
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
            endif
          endif
      end select
!      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
!        rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) + Res
!      endif      
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo

#ifdef CHUAN_CO2
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit

     msrc(:) = source_sink%flow_condition%pressure%dataset%cur_value(:)
     msrc(1) =  msrc(1) / FMWH2O *1D3
     msrc(2) =  msrc(2) / FMWCO2 *1D3
    ! print *,'RT SC source'
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      Res=0D0
      
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS)
           do ieqgas = 1, reaction%ngas
              if(abs(reaction%co2_gas_id) == ieqgas )then
                 icomp = reaction%eqgasspecid(1,ieqgas)
                 iend = local_id*reaction%ncomp
                 istart = iend-reaction%ncomp
                 Res(icomp) = -msrc(2)
                 r_p(istart+icomp) = r_p(istart+icomp) + Res(icomp)
                 print *,'RT SC source', ieqgas,icomp, res(icomp)  
              endif 
           enddo
      end select 
     enddo
     source_sink => source_sink%next
   enddo 
     
#endif
#endif
#if 1  
! Reactions
  if (associated(reaction)) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      Res = 0.d0
      Jup = 0.d0
      call RReaction(Res,Jup,PETSC_FALSE,rt_aux_vars(ghosted_id), &
                     global_aux_vars(ghosted_id), &
                     volume_p(local_id),reaction,option)
      r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)                    

    enddo
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
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

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
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid

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

  ! pass #1 for internal and boundary flux terms  
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
         call SAMRSetCurrentJacobianPatch(J,grid%structured_grid%p_samr_patch)
      endif

      call RTJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  ! pass #2 for everything else
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
         call SAMRSetCurrentJacobianPatch(J,grid%structured_grid%p_samr_patch)
      endif

      call RTJacobianPatch2(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
    
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
      
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Jdn(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc


  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)

  if (option%use_samr) then
  ! Interior Advective Flux Terms -----------------------------------
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

      call TFluxDerivativeAdv(rt_aux_vars(ghosted_id_up), &
                              global_aux_vars(ghosted_id_up), &
                              rt_aux_vars(ghosted_id_dn), &
                              global_aux_vars(ghosted_id_dn), &
                              cur_connection_set%area(iconn), &
                              rt_parameter,option, &
                              patch%internal_velocities(:,iconn),Jup,Jdn)

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

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  ! Boundary Advective Flux Terms -----------------------------------
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

      call TBCFluxDerivativeAdv(boundary_condition%tran_condition%itype, &
                                rt_aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                rt_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                cur_connection_set%area(iconn), &
                                rt_parameter,option, &
                                patch%boundary_velocities(:,sum_connection), &
                                Jdn)
 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  ! Interior Diffusive Flux Terms -----------------------------------
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
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TFluxDerivativeDiff(rt_aux_vars(ghosted_id_up), &
                               global_aux_vars(ghosted_id_up), &
                               porosity_loc_p(ghosted_id_up), &
                               tor_loc_p(ghosted_id_up), &
                               dist_up, &
                               rt_aux_vars(ghosted_id_dn), &
                               global_aux_vars(ghosted_id_dn), &
                               porosity_loc_p(ghosted_id_dn), &
                               tor_loc_p(ghosted_id_dn), &
                               dist_dn, &
                               cur_connection_set%area(iconn), &
                               rt_parameter,option, &
                               patch%internal_velocities(:,iconn),Jup,Jdn)

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

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  ! Boundary Diffusive Flux Terms -----------------------------------
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

      call TBCFluxDerivativeDiff(boundary_condition%tran_condition%itype, &
                                 rt_aux_vars_bc(sum_connection), &
                                 global_aux_vars_bc(sum_connection), &
                                 rt_aux_vars(ghosted_id), &
                                 global_aux_vars(ghosted_id), &
                                 porosity_loc_p(ghosted_id), &
                                 tor_loc_p(ghosted_id), &
                                 cur_connection_set%dist(0,iconn), &
                                 cur_connection_set%area(iconn), &
                                 rt_parameter,option, &
                                 patch%boundary_velocities(:,sum_connection), &
                                 Jdn)
 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  else ! !use_samr

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
      dist_up = distance*fraction_upwind
      dist_dn = distance-dist_up ! should avoid truncation error

      call TFluxDerivative(rt_aux_vars(ghosted_id_up), &
                           global_aux_vars(ghosted_id_up), &
                           porosity_loc_p(ghosted_id_up), &
                           tor_loc_p(ghosted_id_up), &
                           dist_up, &
                           rt_aux_vars(ghosted_id_dn), &
                           global_aux_vars(ghosted_id_dn), &
                           porosity_loc_p(ghosted_id_dn), &
                           tor_loc_p(ghosted_id_dn), &
                           dist_dn, &
                           cur_connection_set%area(iconn), &
                           rt_parameter,option, &
                           patch%internal_velocities(:,iconn),Jup,Jdn)

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

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      call TBCFluxDerivative(boundary_condition%tran_condition%itype, &
                             rt_aux_vars_bc(sum_connection), &
                             global_aux_vars_bc(sum_connection), &
                             rt_aux_vars(ghosted_id), &
                             global_aux_vars(ghosted_id), &
                             porosity_loc_p(ghosted_id), &
                             tor_loc_p(ghosted_id), &
                             cur_connection_set%dist(0,iconn), &
                             cur_connection_set%area(iconn), &
                             rt_parameter,option, &
                             patch%boundary_velocities(:,sum_connection), &
                             Jdn)
 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  endif ! #else AMR_FLUX

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)

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
  
  implicit none

  interface
     subroutine SAMRSetJacobianSrcCoeffsOnPatch(which_pc, p_application, p_patch) 
#include "finclude/petsc.h"

       PetscInt :: which_pc
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSrcCoeffsOnPatch
  end interface
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag  
  type(realization_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), work_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
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
  PetscInt :: iconn
  PetscReal :: qsrc, rdum
  PetscTruth :: volumetric
  
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


  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
    
  if (.not.option%steady_state) then
#if 1  
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      call RTAccumulationDerivative(rt_aux_vars(ghosted_id), &
                                    global_aux_vars(ghosted_id), &
                                    porosity_loc_p(ghosted_id), &
                                    volume_p(local_id),reaction,option,Jup) 
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)                        
    enddo
#endif
  endif
#if 1
  ! Source/Sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)
      if (source_sink%flow_condition%rate%itype == &
          VOLUMETRIC_RATE_SS) then
        volumetric = PETSC_TRUE
      else
        volumetric = PETSC_FALSE
      endif
    endif
      
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      Jup = 0.d0
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          do istart = 1, reaction%ncomp
            Jup(istart,istart) = 1.d-6* &
                                 porosity_loc_p(ghosted_id)* &
                                 global_aux_vars(ghosted_id)%sat(option%liquid_phase)* &
                                 volume_p(local_id)
          enddo
        case(MASS_RATE_SS)
        case default
          if (qsrc < 0) then ! extraction
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              do istart = 1, reaction%ncomp
                Jup(istart,istart) = -qsrc
              enddo
            else
              do istart = 1, reaction%ncomp
                Jup(istart,istart) = -qsrc/global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)
              enddo
            endif
          endif
      end select
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr) 
    enddo                       
    source_sink => source_sink%next
  enddo
#endif
#if 1  
! Reactions
  if (associated(reaction)) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      Res = 0.d0
      Jup = 0.d0
      call RReactionDerivative(Res,Jup,rt_aux_vars(ghosted_id), &
                               global_aux_vars(ghosted_id), &
                               volume_p(local_id),reaction,option)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1, &
                                    Jup,ADD_VALUES,ierr)                        
    enddo
  endif
#endif
 
  if (reaction%use_log_formulation) then
    call GridVecGetArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      iend = ghosted_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) then
          work_loc_p(istart:iend) = 1.d0
        else
          work_loc_p(istart:iend) = rt_aux_vars(ghosted_id)%pri_molal(:)
        endif
      else
        work_loc_p(istart:iend) = rt_aux_vars(ghosted_id)%pri_molal(:)
      endif
    enddo
    call GridVecRestoreArrayF90(grid,field%tran_work_loc, work_loc_p, ierr)
  endif

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_accum, accum_p, ierr)
 
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  
  if (patch%aux%RT%inactive_cells_exist) then
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,rdum,ierr) 
  endif

  if(option%use_samr) then
     tran_pc = 1
     call SAMRSetJacobianSrcCoeffsOnPatch(tran_pc, &
          realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
  endif

end subroutine RTJacobianPatch2

! ************************************************************************** !
!
! RTUpdateAuxVars: Updates the auxilliary variables associated with 
!                  reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVars(realization,update_bcs,update_activity_coefs)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscTruth :: update_bcs
  PetscTruth :: update_activity_coefs
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTUpdateAuxVarsPatch(realization,update_bcs,update_activity_coefs)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine RTUpdateAuxVars

! ************************************************************************** !
!
! RTUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                       reactive transport
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine RTUpdateAuxVarsPatch(realization,update_bcs,compute_activity_coefs)

  use Realization_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  PetscTruth :: update_bcs
  PetscTruth :: compute_activity_coefs
  PetscTruth :: time0
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%reaction%ncomp)
  PetscReal, pointer :: basis_molarity_p(:)
  PetscReal :: weight
  PetscInt, parameter :: iphase = 1
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  reaction => realization%reaction
  
  call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*reaction%ncomp
    istart = iend-reaction%ncomp+1
    
    patch%aux%RT%aux_vars(ghosted_id)%pri_molal = xx_loc_p(istart:iend)
    if (compute_activity_coefs) then
      call RActivityCoefficients(patch%aux%RT%aux_vars(ghosted_id), &
                                 patch%aux%Global%aux_vars(ghosted_id), &
                                 reaction,option)
    endif
    call RTAuxVarCompute(patch%aux%RT%aux_vars(ghosted_id), &
                         patch%aux%Global%aux_vars(ghosted_id), &
                         reaction,option)
  enddo

  if (update_bcs) then
  
    boundary_condition => patch%boundary_conditions%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      basis_molarity_p => boundary_condition%tran_condition% &
        cur_constraint_coupler%aqueous_species%basis_molarity

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif

        select case(boundary_condition%tran_condition%itype)
          case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
            ! since basis_molarity is in molarity, must convert to molality
            ! by dividing by density of water (mol/L -> mol/kg)
            xxbc(1:reaction%ncomp) = basis_molarity_p(1:reaction%ncomp) / &
              patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
          case(DIRICHLET_ZERO_GRADIENT_BC)
!geh            do iphase = 1, option%nphase
              if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                ! same as dirichlet above
                xxbc(1:reaction%ncomp) = basis_molarity_p(1:reaction%ncomp) / &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
              else
                ! same as zero_gradient below
                do idof=1,reaction%ncomp
                  xxbc(idof) = xx_loc_p((ghosted_id-1)*reaction%ncomp+idof)
                enddo
              endif
!geh            enddo
          case(ZERO_GRADIENT_BC)
            do idof=1,reaction%ncomp
              xxbc(idof) = xx_loc_p((ghosted_id-1)*reaction%ncomp+idof)
            enddo
        end select
        ! no need to update boundary fluid density since it is already set
        patch%aux%RT%aux_vars_bc(sum_connection)%pri_molal = xxbc
        if (compute_activity_coefs) then
          call RActivityCoefficients(patch%aux%RT%aux_vars_bc(sum_connection), &
                                     patch%aux%Global%aux_vars_bc(sum_connection), &
                                     reaction,option)
        endif
        call RTAuxVarCompute(patch%aux%RT%aux_vars_bc(sum_connection), &
                             patch%aux%Global%aux_vars_bc(sum_connection), &
                             reaction,option)
      enddo
      boundary_condition => boundary_condition%next
    enddo

    patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE

  endif 
  
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)
  
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
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr

  flag = 0
  grid => patch%grid
  
  n_zero_rows = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + reaction%ncomp
      else
      endif
    enddo
  else
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
        do icomp = 1, reaction%ncomp
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*reaction%ncomp+icomp
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*reaction%ncomp+icomp-1
        enddo
      else
      endif
    enddo
  endif

  patch%aux%RT%zero_rows_local => zero_rows_local
  patch%aux%RT%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%RT%n_zero_rows = n_zero_rows  


  if(.not.(option%use_samr)) then
     call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
          option%mycomm,ierr)
     
     if (flag > 0) patch%aux%RT%inactive_cells_exist = PETSC_TRUE
     
     if (ncount /= n_zero_rows) then
        print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
        stop
     endif
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
! RTGetTecplotHeader: Returns Reactive Transport contribution to 
!                     Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
function RTGetTecplotHeader(realization,icolumn)

  use Realization_module
  use Option_module

  implicit none
  
  character(len=MAXHEADERLENGTH) :: RTGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXHEADERLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  PetscInt :: i
  
  option => realization%option
  reaction => realization%reaction
  
  string = ''
  
  if ((reaction%print_pH) .and. &
      reaction%h_ion_id > 0) then
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-pH"'')') icolumn
    else
      write(string2,'('',"pH"'')') 
    endif
    string = trim(string) // trim(string2)
  endif
  
  if (reaction%print_total_component) then
    do i=1,option%ntrandof
      if (reaction%primary_species_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_tot"'')') icolumn, &
            trim(reaction%primary_species_names(i))
        else
          write(string2,'('',"'',a,''"'')') trim(reaction%primary_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  if (reaction%print_free_ion) then
    do i=1,option%ntrandof
      if (reaction%primary_species_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_free"'')') icolumn, &
            trim(reaction%primary_species_names(i))
        else
          write(string2,'('',"'',a,''"'')') trim(reaction%primary_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo  
  endif
    
  if (reaction%print_act_coefs) then
    do i=1,option%ntrandof
      if (reaction%primary_species_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_gam"'')') icolumn, &
            trim(reaction%primary_species_names(i))
        else
          write(string2,'('',"'',a,''_gam"'')') trim(reaction%primary_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  do i=1,reaction%nkinmnrl
    if (reaction%kinmnrl_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-'',a,''_vf"'')') icolumn, &
          trim(reaction%kinmnrl_names(i))
      else
        write(string2,'('',"'',a,''_vf"'')') trim(reaction%kinmnrl_names(i))    
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,reaction%nkinmnrl
    if (reaction%kinmnrl_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-'',a,''_rt"'')') icolumn, &
          trim(reaction%kinmnrl_names(i))
      else
        write(string2,'('',"'',a,''_rt"'')') trim(reaction%kinmnrl_names(i))    
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,realization%reaction%neqsurfcmplxrxn
    if (reaction%surface_site_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%surface_site_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%surface_site_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,realization%reaction%neqsurfcmplx
    if (reaction%surface_complex_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%surface_complex_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%surface_complex_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo

  if (associated(reaction%kd_print)) then
    do i=1,option%ntrandof
      if (reaction%kd_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_kd"'')') icolumn, &
            trim(reaction%primary_species_names(i))
        else
          write(string2,'('',"'',a,''_kd"'')') trim(reaction%primary_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  if (associated(reaction%total_sorb_print)) then
    do i=1,option%ntrandof
      if (reaction%total_sorb_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_total_sorb"'')') icolumn, &
            trim(reaction%primary_species_names(i))
        else
          write(string2,'('',"'',a,''_total_sorb"'')') trim(reaction%primary_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  RTGetTecplotHeader = string

end function RTGetTecplotHeader

! ************************************************************************** !
!
! RTAuxVarCompute: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 08/28/08
!
! ************************************************************************** !
subroutine RTAuxVarCompute(rt_aux_var,global_aux_var,reaction,option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_aux_var
  type(global_auxvar_type) :: global_aux_var
  
#if 0  
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp
  PetscReal :: dtotal(reaction%ncomp,reaction%ncomp)
  PetscReal :: dtotalsorb(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
#endif

  ! any changes to the below must also be updated in 
  ! Reaction.F90:RReactionDerivative()
  
!already set  rt_aux_var%pri_molal = x

  call RTotal(rt_aux_var,global_aux_var,reaction,option)
  if (reaction%neqsorb > 0) then
    call RTotalSorb(rt_aux_var,global_aux_var,reaction,option)
  endif

#if 0
! numerical check
  Res_orig = 0.d0
  dtotal = 0.d0
  dtotalsorb = 0.d0
  option%iflag = 0 ! be sure not to allocate mass_balance array
  call RTAuxVarInit(rt_auxvar_pert,reaction,option)
  do jcomp = 1, reaction%ncomp
    Res_pert = 0.d0
    call RTAuxVarCopy(rt_auxvar_pert,rt_aux_var,option)
    if (reaction%neqcmplx > 0) then
      aux_var%sec_molal = 0.d0
    endif
    if (reaction%ngas > 0) then
      aux_var%gas_molal = 0.d0
    endif
    if (reaction%neqsurfcmplxrxn > 0) then
      rt_auxvar_pert%eqsurfcmplx_freesite_conc = 1.d-9
      rt_auxvar_pert%eqsurfcmplx_conc = 0.d0
    endif
    if (reaction%neqionxrxn > 0) then
      rt_aux_var%eqionx_ref_cation_sorbed_conc = 1.d-9
    endif
    pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
    rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
    
    call RTotal(rt_auxvar_pert,global_aux_var,reaction,option)
    dtotal(:,jcomp) = (rt_auxvar_pert%total(:,1) - rt_aux_var%total(:,1))/pert
    if (reaction%neqsorb > 0) then
      call RTotalSorb(rt_auxvar_pert,global_aux_var,reaction,option)
      if (reaction%kinmr_nrate <= 0) &
        dtotalsorb(:,jcomp) = (rt_auxvar_pert%total_sorb_eq(:) - &
                               rt_aux_var%total_sorb_eq(:))/pert
    endif
  enddo
  do icomp = 1, reaction%ncomp
    do jcomp = 1, reaction%ncomp
      if (dabs(dtotal(icomp,jcomp)) < 1.d-16) dtotal(icomp,jcomp) = 0.d0
      if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
        if (dabs(dtotalsorb(icomp,jcomp)) < 1.d-16) dtotalsorb(icomp,jcomp) = 0.d0
      endif
    enddo
  enddo
  rt_aux_var%dtotal(:,:,1) = dtotal
  if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) &
    rt_aux_var%dtotal_sorb = dtotalsorb
  call RTAuxVarDestroy(rt_auxvar_pert)
#endif
  
end subroutine RTAuxVarCompute

! ************************************************************************** !
!
! RTDestroy: Deallocates variables associated with Reactive Transport
! author: Glenn Hammond
! date: 02/03/09
!
! ************************************************************************** !
subroutine RTDestroy(realization)

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
      call RTDestroyPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  
end subroutine RTDestroyPatch

end module Reactive_Transport_module
