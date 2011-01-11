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
            RTGetTecplotHeader, &
            RTUpdateAuxVars, &
            RTComputeMassBalance, &
            RTDestroy, &
            RTUpdateTransportCoefs, &
            RTUpdateRHSCoefs, &
            RTCalculateRHS_t0, &
            RTCalculateRHS_t1, &
            RTCalculateTransportMatrix, &
            RTReact, &
            RTTransportResidual, &
            RTTransportMatVec, &
            RTCheckUpdate, &
            RTJumpStartKineticSorption, &
            RTCheckpointKineticSorption
  
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
  patch%aux%RT%rt_parameter%ncomp = reaction%ncomp
  patch%aux%RT%rt_parameter%naqcomp = reaction%naqcomp
  patch%aux%RT%rt_parameter%nimcomp = 0
  patch%aux%RT%rt_parameter%offset_aq = reaction%offset_aq
  if (reaction%ncollcomp > 0) then
    patch%aux%RT%rt_parameter%ncoll = reaction%ncoll
    patch%aux%RT%rt_parameter%offset_coll = reaction%offset_coll
    patch%aux%RT%rt_parameter%ncollcomp = reaction%ncollcomp
    patch%aux%RT%rt_parameter%offset_collcomp = reaction%offset_collcomp
    allocate(patch%aux%RT%rt_parameter%pri_spec_to_coll_spec(reaction%naqcomp))
    patch%aux%RT%rt_parameter%pri_spec_to_coll_spec = &
      reaction%pri_spec_to_coll_spec
    allocate(patch%aux%RT%rt_parameter%coll_spec_to_pri_spec(reaction%ncollcomp))
    patch%aux%RT%rt_parameter%coll_spec_to_pri_spec = &
      reaction%coll_spec_to_pri_spec
  endif
    
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
  
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%RT%aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RTAuxVarInit(patch%aux%RT%aux_vars_bc(iconn),reaction,option)
    enddo
#ifdef DASVYAT
!	write(*,*) 'sum_connection',sum_connection
!    do iconn = 1, sum_connection
!      write(*,*) "RTAuxVarInit ",patch%aux%RT%aux_vars_bc(iconn)%total(1,1)
!    enddo
!	read(*,*)
#endif
  endif
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
  PetscBool :: changed
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
    if (patch%imat(ghosted_id) <= 0) cycle
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
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy,ierr)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%tran_xx, &
                                   realization%field%tran_xx_loc,NTRANDOF)
  
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
  PetscInt :: ghosted_id, local_id, imnrl, iaqspec, ncomp, icomp
  PetscInt :: k, irate, irxn, icplx, ncplx
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  rt_aux_vars => patch%aux%RT%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

  ! update:                             cells      bcs         act. coefs.
  call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

  if (.not.option%init_stage) then
    ! update mineral volume fractions
    if (reaction%nkinmnrl > 0) then
    
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle
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

#ifdef CHUAN_CO2
          if(option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE)then
            ncomp = reaction%kinmnrlspecid(0,imnrl)
            do iaqspec=1, ncomp  
              icomp = reaction%kinmnrlspecid(iaqspec,imnrl)
              if(icomp == realization%reaction%species_idx%co2_aq_id) then
                global_aux_vars(ghosted_id)%reaction_rate(2) &
                  = global_aux_vars(ghosted_id)%reaction_rate(2)& 
                  + rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* option%tran_dt&
                  * reaction%mnrlstoich(icomp,imnrl)/option%flow_dt
              else if(icomp == reaction%species_idx%h2o_aq_id)then
                global_aux_vars(ghosted_id)%reaction_rate(1) &
                  = global_aux_vars(ghosted_id)%reaction_rate(1)& 
                  + rt_aux_vars(ghosted_id)%mnrl_rate(imnrl)* option%tran_dt&
                  * reaction%mnrlstoich(icomp,imnrl)/option%flow_dt
              endif
            enddo 
          endif   
#endif

        enddo
      enddo
    endif
  

    ! update multirate sorption concentrations 
  ! WARNING: below assumes site concentration multiplicative factor
    if (reaction%kinmr_nrate > 0) then 
      do ghosted_id = 1, grid%ngmax 
        if (patch%imat(ghosted_id) <= 0) cycle
        do irate = 1, reaction%kinmr_nrate 
          kdt = reaction%kinmr_rate(irate) * option%tran_dt 
          one_plus_kdt = 1.d0 + kdt 
          k_over_one_plus_kdt = reaction%kinmr_rate(irate)/one_plus_kdt 
          rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate) = & 
            (rt_aux_vars(ghosted_id)%kinmr_total_sorb(:,irate) + & 
            kdt * reaction%kinmr_frac(irate) * &
            rt_aux_vars(ghosted_id)%total_sorb_eq)/one_plus_kdt
        enddo 
      enddo 
    endif

    ! update kinetic sorption concentrations
    if (reaction%nkinsrfcplxrxn > 0) then
      do ghosted_id = 1, grid%ngmax 
        if (patch%imat(ghosted_id) <= 0) cycle
        do irxn = 1, reaction%nkinsrfcplxrxn
          ncplx = reaction%kinsrfcplx_rxn_to_complex(0,irxn)
          do k = 1, ncplx ! ncplx in rxn
            icplx = reaction%kinsrfcplx_rxn_to_complex(k,irxn)
            rt_aux_vars(ghosted_id)%kinsrfcplx_conc(icplx) = &
              rt_aux_vars(ghosted_id)%kinsrfcplx_conc_kp1(icplx)
          enddo
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
  PetscInt :: dof_offset, istart, iendaq, iendall
  PetscInt :: istartcoll, iendcoll
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
  call GridVecGetArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  call GridVecGetArrayF90(grid,field%tran_accum, accum_p, ierr)

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
      istartcoll = dof_offset + reaction%offset_coll + 1
      iendcoll = dof_offset + reaction%offset_coll + reaction%ncoll
      rt_aux_vars(ghosted_id)%colloid%conc_mob = xx_p(istartcoll:iendcoll)* &
        global_aux_vars(ghosted_id)%den_kg(1)*1.d-3
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
                        reaction,option,accum_p(istart:iendall)) 
    if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
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
      call RTUpdateTransportCoefsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
                      tor_loc_p(ghosted_id_up),dist_up, &
                      global_aux_vars(ghosted_id_dn), &
                      porosity_loc_p(ghosted_id_dn), &
                      tor_loc_p(ghosted_id_dn),dist_dn, &
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
        local_id = cur_connection_set%id_dn(id)
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
      call RTUpdateRHSCoefsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
      call RTCalculateRHS_t0Patch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
      call RTCalculateRHS_t1Patch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aq + 1
      iendaq = reaction%offset_aq + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_coll + 1
        iendcoll = reaction%offset_coll + reaction%ncoll
      endif
      
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          ! units should be mol/sec
          Res(istartaq:iendaq) = -1.d-6* &
                porosity_loc_p(ghosted_id)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
!                (source_sink%tran_condition%cur_constraint_coupler% &
!                 rt_auxvar%total(:,iphase) - rt_aux_vars(ghosted_id)%total(:,iphase))* &
! keeping only the fixed portion for RHS
                (source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%total(:,iphase))* & !- rt_aux_vars(ghosted_id)%total(:,iphase))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
          if (reaction%ncoll > 0) then
            Res(istartcoll:iendcoll)  =  -1.d-6* &
              !in this case, conc_mob is in molality 
              ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)* 
              !         (m^3 bulk)*(kg water/m^3 water)/(sec)*
              !         (mol colloid /kg water) = mol colloid/sec
                porosity_loc_p(ghosted_id)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
!                (source_sink%tran_condition%cur_constraint_coupler% &
!                 rt_auxvar%colloid%conc_mob(:) - rt_aux_vars(ghosted_id)%colloid%conc_mob(:))* &
! keeping only the fixed portion for RHS
                (source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%colloid%conc_mob(:))* & ! - rt_aux_vars(ghosted_id)%colloid%conc_mob(:))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
          endif
        case(MASS_RATE_SS)
          Res(istartaq:iendaq) = -source_sink%tran_condition% &
                 cur_constraint_coupler%rt_auxvar%total(:,iphase) ! actually moles/sec
          if (reaction%ncoll > 0) then
            option%io_buffer = 'Need to implement MASS_RATE_SS source/sink ' // &
                               'term correctly for colloids'
            call printErrMsg(option)
            Res(istartcoll:iendcoll) = -source_sink%tran_condition% &
                 cur_constraint_coupler%rt_auxvar%colloid%conc_mob(:) ! actually moles/sec
          endif          
        case default
          if (qsrc > 0) then ! injection
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              Res(istartaq:iendaq) = -qsrc* & ! m^3 water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%total(:,iphase)*1000.d0
              if (reaction%ncoll > 0) then
                Res(istartcoll:iendcoll) = -qsrc* & ! m^3 water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%colloid%conc_mob(:)*1000.d0
              endif                     
            else ! mass
              Res(istartaq:iendaq) = -qsrc* & ! kg water / sec
                     source_sink%tran_condition%cur_constraint_coupler% &
                     rt_auxvar%total(:,iphase)/ &
                     global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                     1000.d0
              if (reaction%ncoll > 0) then  ! needs to be moles/sec
                Res(istartcoll:iendcoll) = -qsrc* & ! kg water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%colloid%conc_mob(:)/ &
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0
              endif                     
            endif
          else ! extraction
            ! all extraction, which is a function of the cell concentration is
            ! handled in the transport matrix, not the RHS
!            if (volumetric) then ! qsrc is volumetric; must be converted to mass
!              Res(istartaq:iendaq) = -qsrc*rt_aux_vars(ghosted_id)%total(:,iphase)*1000.d0
!              if (reaction%ncoll > 0) then
!                Res(istartcoll:iendcoll) = -qsrc* & ! m^3 water / sec
!                    rt_aux_vars(ghosted_id)%colloid%conc_mob(:)*1000.d0
!              endif               
!            else
!              Res(istartaq:iendaq) = -qsrc* &
!                    rt_aux_vars(ghosted_id)%total(:,iphase)/ &
!                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
!                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
!              if (reaction%ncoll > 0) then  ! needs to be moles/sec
!                Res(istartcoll:iendcoll) = -qsrc* & ! kg water / sec
!                    rt_aux_vars(ghosted_id)%colloid%conc_mob(:)/ & 
!                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
!                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
!              endif                     
!            endif
          endif
      end select
!      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
!        rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) + Res
!      endif
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      rhs_p(istartall:iendall) = rhs_p(istartall:iendall) + Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo

#ifdef CHUAN_CO2
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      source_sink => patch%source_sinks%first 
      do 
        if (.not.associated(source_sink)) exit

        msrc(:) = source_sink%flow_condition%pressure%dataset%cur_value(:)
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
                if(abs(reaction%species_idx%co2_gas_id) == ieqgas) then
                  icomp = reaction%eqgasspecid(1,ieqgas)
                  iendall = local_id*reaction%ncomp
                  istartall = iendall-reaction%ncomp
                  Res(icomp) = -msrc(2)
                  rhs_p(istartall+icomp) = rhs_p(istartall+icomp) + Res(icomp)
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
  use Level_module
  use Patch_module
  use Grid_module

  implicit none

  interface
     subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "finclude/petscsysdef.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
       
       Mat :: mat
       PetscFortranAddr :: patch
     end subroutine SAMRSetCurrentJacobianPatch
  end interface
      
  type(realization_type) :: realization
  Mat :: T
  
  type(grid_type), pointer :: grid
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option 
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  option => realization%option
 
  call MatZeroEntries(T,ierr)
  
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
      if(option%use_samr) then
         call SAMRSetCurrentJacobianPatch(T, grid%structured_grid%p_samr_patch)
      endif

      call RTCalculateTranMatrixPatch1(realization,T)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
      if(option%use_samr) then
         call SAMRSetCurrentJacobianPatch(T, grid%structured_grid%p_samr_patch)
      endif

      call RTCalculateTranMatrixPatch2(realization,T)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
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
  PetscReal :: qsrc
  PetscBool :: volumetric  
  PetscErrorCode :: ierr
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
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
                     coef_up,coef_dn)

 !     coef_dn = coef_dn*global_aux_vars(ghosted_id)%den_kg*1.d-3

      !Jup not needed 
      coef_dn = -coef_dn
      call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr)
    
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! need to add source/sink
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

      if (patch%imat(ghosted_id) <= 0) cycle
      
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          ! units should be mol/sec
          coef_dn(1) = -1.d-6* &
                porosity_loc_p(ghosted_id)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
! keeping only the non-fixed portion for RHS
                !(source_sink%tran_condition%cur_constraint_coupler% &
                !rt_auxvar%total(:,iphase) 
                -1.d0 * & !rt_aux_vars(ghosted_id)%total(:,iphase))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
        case(MASS_RATE_SS)
          ! since specified mass rate is fixed, include only on RHS
        case default
          if (qsrc > 0) then ! injection
            ! all injection, which is fixed (i.e. not a function of the 
            ! cell concentration is handled on the RHS, not in the transport 
            ! matrix
          else ! extraction
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              coef_dn(1) = -qsrc*1000.d0 !rt_aux_vars(ghosted_id)%total(:,iphase)*1000.d0
            else
              coef_dn(1) = -qsrc* &
                    1.d0/ & !rt_aux_vars(ghosted_id)%total(:,iphase)/ &
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
            endif
          endif
      end select
      call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr)

    enddo
    source_sink => source_sink%next
  enddo

  ! All CO2 source/sinks are handled on the RHS for now
#endif

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
 
  interface

     subroutine SAMRSetJacobianSrcCoeffsOnPatch(which_pc, p_application, p_patch) 
#include "finclude/petscsysdef.h"

       PetscInt :: which_pc
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSrcCoeffsOnPatch
  end interface
 
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
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(1), coef_dn(1)
  PetscErrorCode :: ierr
  PetscInt :: flow_pc
    
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
                        
  ! need to add source/sink

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

  if(option%use_samr) then
     flow_pc = 1
     call SAMRSetJacobianSrcCoeffsOnPatch(flow_pc, &
          realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
  endif
  
end subroutine RTCalculateTranMatrixPatch2

! ************************************************************************** !
!
! RTReactPatch: Calculate reaction
! author: Glenn Hammond
! date: 05/03/10
!
! ************************************************************************** !
subroutine RTReact(realization)

  use Realization_module
  use Field_module
  use Discretization_module    
  use Level_module
  use Patch_module
  use Option_module
  use Logging_module

  implicit none
      
  interface
     subroutine SAMRCoarsenVector(p_application, vec)
       implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: p_application
       Vec :: vec
      end subroutine SAMRCoarsenVector
  end interface 

  type(realization_type) :: realization
  type(discretization_type), pointer :: discretization
  type(field_type), pointer ::field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
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

  if(option%use_samr) then
    call SAMRCoarsenVector(discretization%amrgrid%p_application, field%tran_xx)
  endif
      
#ifdef OS_STATISTICS
  call_count = 0
  sum_newton_iterations = 0
  max_newton_iterations_in_a_cell = -99999999
  max_newton_iterations_on_a_core = -99999999
  min_newton_iterations_on_a_core = 99999999
#endif  

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
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

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
     
  implicit none
  
  type(realization_type) :: realization
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscInt :: iphase
  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: mask_p(:)
#ifdef CHUNK
  PetscInt :: num_iterations(realization%option%chunk_size)
  PetscInt :: local_start
  PetscInt :: local_end
  PetscInt :: ghosted_start
  PetscInt :: ghosted_end
  PetscInt :: chunk_size_save
  PetscInt :: ichunk
#else
  PetscInt :: num_iterations
#endif
#ifdef OS_STATISTICS
  PetscInt :: sum_iterations
  PetscInt :: max_iterations
  PetscInt :: icount
#endif
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  rt_aux_vars => patch%aux%RT%aux_vars
  grid => patch%grid
  reaction => realization%reaction

  ! need up update aux vars based on current density/saturation,
  ! but NOT activity coefficients
  call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

  ! Get vectors
  call GridVecGetArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)

  if(option%use_samr) then
      call GridVecGetMaskArrayCellF90(grid, field%porosity_loc, mask_p, ierr)
  endif
      
  iphase = 1
#ifdef OS_STATISTICS
  sum_iterations = 0
  max_iterations = 0
  icount = 0
#endif

#ifdef CHUNK
  local_start = 1
  chunk_size_save = option%chunk_size
  do
    ! change chunk size if array is not evenly divisble by chunk_size
    if (local_start + option%chunk_size - 1 > grid%nlmax) then
      option%chunk_size = grid%nlmax - local_start + 1
    endif

    local_end = local_start+option%chunk_size-1
  
    ghosted_start = grid%nL2G(local_start)
    ghosted_end = grid%nL2G(local_end)

    ! We actually need to check all cells for inactives
    if (patch%imat(ghosted_start) <= 0) cycle

    iend = local_end*reaction%naqcomp
    istart = iend-reaction%naqcomp*option%chunk_size+1

    ! tran_xx_p passes in total component concentrations
    !       and returns free ion concentrations
    call RReactChunk(rt_aux_vars(ghosted_start:ghosted_end), &
                global_aux_vars(ghosted_start:ghosted_end), &
                tran_xx_p(istart:iend),volume_p(local_start:local_end), &
                porosity_loc_p(ghosted_start:ghosted_end), &
                num_iterations,reaction,option)
    ! set primary dependent var back to free-ion molality
    ! NOW THE BELOW IS PERFORMED WITHIN RReactChunk()
    !tran_xx_p(istart:iend) = rt_aux_vars(ghosted_id)%pri_molal
#ifdef OS_STATISTICS
    do ichunk = 1, option%chunk_size
      if (num_iterations(ichunk) > max_iterations) then
        max_iterations = num_iterations(ichunk)
      endif
      sum_iterations = sum_iterations + num_iterations(ichunk)
      icount = icount + 1
    enddo
#endif
 
   if (local_end >= grid%nlmax) then
     option%chunk_size = chunk_size_save
     exit
   else
     local_start = local_start + option%chunk_size
   endif
 
  enddo
#else
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if(option%use_samr .and. (mask_p(local_id)<=0)) cycle
      
    iend = local_id*reaction%naqcomp
    istart = iend-reaction%naqcomp+1
    call RReact(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                tran_xx_p(istart:iend),volume_p(local_id), &
                porosity_loc_p(ghosted_id), &
                num_iterations,reaction,option)
    ! set primary dependent var back to free-ion molality
    tran_xx_p(istart:iend) = rt_aux_vars(ghosted_id)%pri_molal
#ifdef OS_STATISTICS
    if (num_iterations > max_iterations) then
      max_iterations = num_iterations
    endif
    sum_iterations = sum_iterations + num_iterations
    icount = icount + 1
#endif
  enddo
#endif  
  
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
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  
  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)

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

end subroutine RTComputeBCMassBalanceOSPatch

! ************************************************************************** !
!
! RTTransportResidual: Calculates the transport residual equation for a single 
!                      chemical component (total component concentration)
! author: Glenn Hammond
! date: 05/24/10
!
! ************************************************************************** !
subroutine RTTransportResidual(realization,solution_loc,residual,idof)

  use Realization_module
  use Field_module
  use Level_module
  use Patch_module
  use Discretization_module
  use Option_module
      
  interface
  subroutine samrpetscobjectstateincrease(vec)
  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  Vec :: vec
  end subroutine samrpetscobjectstateincrease
  end interface

  type(realization_type) :: realization
  Vec :: solution_loc
  Vec :: residual
  PetscInt :: idof
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  type(option_type), pointer :: option
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTTransportResidualPatch1(realization,solution_loc,residual,idof)
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
        call RTTransportResidualFluxContribPatch(residual,realization,ierr)
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
      call RTTransportResidualPatch2(realization,solution_loc,residual,idof)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(discretization%itype==AMR_GRID) then
    call samrpetscobjectstateincrease(residual)
  endif
      
end subroutine RTTransportResidual

! ************************************************************************** !
!
! RTTransportResidualPatch: Calculates the transport residual equation for  
!                           a single chemical component (total component 
!                           concentration)
! author: Glenn Hammond
! date: 05/24/10
!
! ************************************************************************** !
subroutine RTTransportResidualPatch1(realization,solution_loc,residual,idof)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none

  interface
     PetscInt function samr_patch_at_bc(p_patch, axis, dim)
     implicit none
     
#include "finclude/petscsysdef.h"
     
     PetscFortranAddr :: p_patch
     PetscInt :: axis,dim
   end function samr_patch_at_bc
  end interface

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  
  type(realization_type) :: realization
  Vec :: solution_loc
  Vec :: residual
  PetscInt :: idof
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: solution_loc_p(:)
  PetscReal, pointer :: residual_p(:)
  PetscReal, pointer :: rhs_coef_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  PetscReal :: res
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(1), coef_dn(1)
  PetscErrorCode :: ierr

  ! samr specific
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecGetArrayF90(grid,residual, residual_p, ierr)  
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  

  if (option%use_samr) then
     do axis=0,2  
        call GridVecGetArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
     enddo
  endif

  iphase = 1
  
  residual_p = 0.d0

  if (option%use_samr) then
    nlx = grid%structured_grid%nlx  
    nly = grid%structured_grid%nly  
    nlz = grid%structured_grid%nlz 

    ngx = grid%structured_grid%ngx   
    ngxy = grid%structured_grid%ngxy

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, ZERO_INTEGER, &
                        ZERO_INTEGER)==1) nlx = nlx-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, ZERO_INTEGER, &
                        ONE_INTEGER)==1) nlx = nlx-1
    
    max_x_conn = (nlx+1)*nly*nlz
    ! reinitialize nlx
    nlx = grid%structured_grid%nlx  

    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, ONE_INTEGER, &
                        ZERO_INTEGER)==1) nly = nly-1
    if(samr_patch_at_bc(grid%structured_grid%p_samr_patch, ONE_INTEGER, &
                        ONE_INTEGER)==1) nly = nly-1
    
    max_y_conn = max_x_conn + nlx*(nly+1)*nlz

    ! reinitialize nly
    nly = grid%structured_grid%nly  
  endif

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
                     coef_up,coef_dn)
      res = coef_up(iphase)*solution_loc_p(ghosted_id_up) + &
            coef_dn(iphase)*solution_loc_p(ghosted_id_dn)

      
      if (option%use_samr) then
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
        fluxes(direction)%flux_p(flux_id) = res
      endif

      if(.not.option%use_samr) then
        residual_p(local_id_up) = residual_p(local_id_up) + res
        residual_p(local_id_dn) = residual_p(local_id_dn) - res
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
                     coef_up,coef_dn)

      ! leave off the boundary contribution since it is already inclued in the
      ! rhs value
      res = coef_up(iphase)*rt_aux_vars_bc(sum_connection)%total(idof,iphase) + &
            coef_dn(iphase)*solution_loc_p(ghosted_id)
!geh - the below assumes that the boundary contribution was provided in the
!      rhs vector above.
!geh      res = coef_dn(iphase)*solution_loc_p(ghosted_id)

      if (option%use_samr) then
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
              fluxes(direction)%flux_p(flux_id) = res
           case(EAST_FACE)
              ghosted_id = ghosted_id+1
              flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)
              fluxes(direction)%flux_p(flux_id) = -res
           case(SOUTH_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = res
           case(NORTH_FACE)
              ghosted_id = ghosted_id+ngx
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -res
           case(BOTTOM_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = res
           case(TOP_FACE)
              ghosted_id = ghosted_id+ngxy
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -res
         end select
      else
         residual_p(local_id)= residual_p(local_id) - res
      endif
    
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  
  call GridVecRestoreArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,residual, residual_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

end subroutine RTTransportResidualPatch1

! ************************************************************************** !
!
! RTTransportResidualPatch: Calculates the transport residual equation for  
!                           a single chemical component (total component 
!                           concentration)
! author: Glenn Hammond
! date: 05/24/10
!
! ************************************************************************** !
subroutine RTTransportResidualPatch2(realization,solution_loc,residual,idof)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  Vec :: solution_loc
  Vec :: residual
  PetscInt :: idof
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: solution_loc_p(:)
  PetscReal, pointer :: residual_p(:)
  PetscReal, pointer :: rhs_coef_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  PetscReal :: res
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(1), coef_dn(1)
  PetscErrorCode :: ierr
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecGetArrayF90(grid,residual, residual_p, ierr)  
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  
  
  iphase = 1

  ! need to add source/sink

  ! Accumulation term
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    coef = porosity_loc_p(ghosted_id)* &
           global_aux_vars(ghosted_id)%sat(iphase)* &
           1000.d0* &
           volume_p(local_id)/option%tran_dt
!geh need to separate out accumulation from boundary fluxes
!geh    residual_p(local_id) = coef*solution_loc_p(ghosted_id)- &
!geh ! RHS will not change, but is moved to lhs (thus, the -1.d0*)
!geh                           rhs_p((local_id-1)*option%ntrandof+idof)

!geh    residual_p(local_id) = coef*solution_loc_p(ghosted_id)- &

    residual_p(local_id) = residual_p(local_id) + &
                           coef*solution_loc_p(ghosted_id)- &
                           rhs_coef_p(local_id)* &
                           rt_aux_vars(ghosted_id)%total(idof,iphase)

  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  
  call GridVecRestoreArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,residual, residual_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

end subroutine RTTransportResidualPatch2

! ************************************************************************** !
!
! RTTransportMatVecPatch2: Calculates the transport residual equation for  
!                           a single chemical component (total component 
!                           concentration)
! author: Bobby Philip
! date: 06/11/2010
!
! ************************************************************************** !
subroutine RTTransportMatVecPatch2(realization,solution_loc,residual,idof)

  use Realization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_type) :: realization
  Vec :: solution_loc
  Vec :: residual
  PetscInt :: idof
  
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars_bc(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: solution_loc_p(:)
  PetscReal, pointer :: residual_p(:)
  PetscReal, pointer :: rhs_coef_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  PetscReal :: res
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(1), coef_dn(1)
  PetscErrorCode :: ierr
    
  ! Get vectors
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_aux_vars => patch%aux%Global%aux_vars
  rt_aux_vars => patch%aux%RT%aux_vars
  rt_aux_vars_bc => patch%aux%RT%aux_vars_bc
  grid => patch%grid

  ! Get vectors
  call GridVecGetArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecGetArrayF90(grid,residual, residual_p, ierr)  
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  
  
  iphase = 1

  ! need to add source/sink

  ! Accumulation term
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    coef = porosity_loc_p(ghosted_id)* &
           global_aux_vars(ghosted_id)%sat(iphase)* &
           1000.d0* &
           volume_p(local_id)/option%tran_dt
!geh need to separate out accumulation from boundary fluxes
!geh    residual_p(local_id) = coef*solution_loc_p(ghosted_id)- &
!geh ! RHS will not change, but is moved to lhs (thus, the -1.d0*)
!geh                           rhs_p((local_id-1)*option%ntrandof+idof)

!geh    residual_p(local_id) = coef*solution_loc_p(ghosted_id)- &

    residual_p(local_id) = residual_p(local_id) + &
                           coef*solution_loc_p(ghosted_id)

  enddo

  ! Restore vectors
  call GridVecRestoreArrayF90(grid,field%tran_rhs_coef,rhs_coef_p,ierr)  
  call GridVecRestoreArrayF90(grid,solution_loc, solution_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,residual, residual_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)

end subroutine RTTransportMatVecPatch2

subroutine RTTransportMatVec(mat, x, y)

!  use Simulation_module
!  use Timestepper_module
!  use Solver_module
  use Realization_module
  use Discretization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Debug_module
  use ISO_C_BINDING
      
  implicit none

#ifndef PC_BUG
  interface
    subroutine SAMRGetRealization(p_application, realization) 
      use Realization_module
#include "finclude/petscsys.h"
      PetscFortranAddr :: p_application
      type(realization_type), pointer :: realization
      end subroutine SAMRGetRealization

     subroutine SAMRGetPetscTransportMatrix(p_application, transportMat) 
      use Realization_module
#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
      
      PetscFortranAddr :: p_application
      Mat :: transportMat
      end subroutine SAMRGetPetscTransportMatrix      

  end interface
#endif

  Mat, intent(in) :: mat    
  Vec, intent(in) :: x
  Vec, intent(out) :: y

!  type(simulation_type),pointer :: simulation
  type(realization_type),pointer :: realization
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(c_ptr) :: realization_cptr
  PetscInt :: idof
  PetscReal :: alpha    
  PetscErrorCode :: ierr
  PetscFortranAddr :: p_application
  PetscReal :: diff
  Mat :: vmat    
      
  call MatShellGetContext(mat, p_application, ierr)
#ifndef PC_BUG  
  call SAMRGetRealization(p_application, realization)
#endif
!  call SAMRGetPetscTransportMatrix(p_application, vmat)

  field => realization%field
  discretization => realization%discretization
  option => realization%option
  idof = option%rt_idof
      
  ! solution is stored in x vector, but we need it in ghosted
  ! form for the residual calculation
  call DiscretizationGlobalToLocal(discretization,x, &
                                       field%work_loc,ONEDOF)

   ! temporary for samr testing                                     
  call DiscretizationGlobalToLocal(discretization,x, &
                                       field%work_samr_loc,ONEDOF)

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTTransportResidualPatch1(realization,field%work_loc,y,idof)
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
        call RTTransportResidualFluxContribPatch(y,realization,ierr)
        cur_patch => cur_patch%next
      enddo
      cur_level => cur_level%next
    enddo
      endif

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RTTransportMatVecPatch2(realization,field%work_loc,y,idof)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

!  alpha=-1.0    
!  call VecScale(y, alpha, ierr)

!  call MatMult(vmat, x, field%work_samr, ierr)

!  call VecAXPY(field%work_samr, alpha, y, ierr)

!    call VecNorm(field%work_samr, NORM_2, diff, ierr)
      
end subroutine RTTransportMatVec
      
! ************************************************************************** !
!
! RichardsResidualFluxContribsPatch: should be called only for SAMR
! author: Bobby Philip
! date: 02/17/09
!
! ************************************************************************** !
subroutine RTTransportResidualFluxContribPatch(r,realization,ierr)
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
  PetscInt :: axis, nlx, nly, nlz
  PetscInt :: iconn, i, j, k
  PetscInt :: xup_id, xdn_id, yup_id, ydn_id, zup_id, zdn_id

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)

  do axis=0,2  
     call GridVecGetArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
  enddo

  nlx = grid%structured_grid%nlx  
  nly = grid%structured_grid%nly  
  nlz = grid%structured_grid%nlz 
  
  iconn=0
  do k=1,nlz
     do j=1,nly
        do i=1,nlx
           iconn=iconn+1
           xup_id = ((k-1)*nly+j-1)*(nlx+1)+i
           xdn_id = xup_id+1
           yup_id = ((k-1)*(nly+1)+(j-1))*nlx+i
           ydn_id = yup_id+nlx
           zup_id = ((k-1)*nly+(j-1))*nlx+i
           zdn_id = zup_id+nlx*nly

           r_p(iconn) = r_p(iconn)+fluxes(0)%flux_p(xdn_id)-fluxes(0)%flux_p(xup_id) &
                                  +fluxes(1)%flux_p(ydn_id)-fluxes(1)%flux_p(yup_id) &
                                  +fluxes(2)%flux_p(zdn_id)-fluxes(2)%flux_p(zup_id)

        enddo
     enddo
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)

end subroutine RTTransportResidualFluxContribPatch

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
  use Level_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Logging_module

  implicit none
  
  interface
  subroutine samrpetscobjectstateincrease(vec)
  implicit none
#include "finclude/petscsys.h"
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
  
  call PetscLogEventBegin(logging%event_rt_residual,ierr)

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
  
  call PetscLogEventEnd(logging%event_rt_residual,ierr)

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
     
#include "finclude/petscsys.h"
     
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
  
  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)

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
  
  if (.not.patch%aux%RT%aux_vars_up_to_date) then
    if (reaction%act_coef_update_frequency == ACT_COEF_FREQUENCY_NEWTON_ITER) then
      ! update:                             cells      bcs        act. coefs.
      call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
    else 
      ! update:                             cells      bcs        act. coefs.
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

  if (option%use_samr) then
    ! this section should be broken due to changes in the transport data structures
!!!!!!!!!!!!!!!!#ifndef REVISED_TRANSPORT
#if 0
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
          
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle

        call TFluxAdv(rt_aux_vars(ghosted_id_up),global_aux_vars(ghosted_id_up), &
               rt_aux_vars(ghosted_id_dn),global_aux_vars(ghosted_id_dn), &
               cur_connection_set%area(iconn),rt_parameter,option, &
               patch%internal_velocities(:,sum_connection),Res)
          
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
          
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
          
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
               dist_dn, &
               cur_connection_set%area(iconn),rt_parameter,option, &
               patch%internal_velocities(:,sum_connection),Res)
          
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
          
        if (patch%imat(ghosted_id) <= 0) cycle
          
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
          
        if (patch%imat(ghosted_id) <= 0) cycle
          
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
          
        direction = (boundary_condition%region%faces(iconn)-1)/2
          
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
!!!!!!!!!!!!#endif ! #ifndef REVISED_TRANSPORT
#endif
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

        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle

        ! TFluxCoef will eventually be moved to another routine where it should be
        ! called only once per flux interface at the beginning of a transport
        ! time step.
        call TFluxCoef(option,cur_connection_set%area(iconn), &
                       patch%internal_velocities(:,sum_connection), &
                       patch%internal_tran_coefs(:,sum_connection), &
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

        ! TFluxCoef accomplishes the same as what TBCCoef would
        call TFluxCoef(option,cur_connection_set%area(iconn), &
                       patch%boundary_velocities(:,sum_connection), &
                       patch%boundary_tran_coefs(:,sum_connection), &
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
      
      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif

#ifdef DASVYAT
!  do iconn = 1, grid%nlmax
!	  write(*,*) "r_p", r_p(iconn)
!  end do
!  stop
#endif  
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
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:), rt_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscReal :: qsrc, molality
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscBool :: volumetric
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
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp

      call RTAccumulation(rt_aux_vars(ghosted_id), &
                          global_aux_vars(ghosted_id), &
                          porosity_loc_p(ghosted_id), &
                          volume_p(local_id), &
                          reaction,option,Res)
      if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
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

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aq + 1
      iendaq = reaction%offset_aq + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_coll + 1
        iendcoll = reaction%offset_coll + reaction%ncoll
      endif
      
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          ! units should be mol/sec
          Res(istartaq:iendaq) = -1.d-6* &
                porosity_loc_p(ghosted_id)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
                (source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%total(:,iphase) - rt_aux_vars(ghosted_id)%total(:,iphase))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
          if (reaction%ncoll > 0) then
            Res(istartcoll:iendcoll)  =  -1.d-6* &
              !in this case, conc_mob is in molality 
              ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)* 
              !         (m^3 bulk)*(kg water/m^3 water)/(sec)*
              !         (mol colloid /kg water) = mol colloid/sec
                porosity_loc_p(ghosted_id)* &
                volume_p(local_id)* & ! convert m^3 water -> L water
                (source_sink%tran_condition%cur_constraint_coupler% &
                 rt_auxvar%colloid%conc_mob(:) - rt_aux_vars(ghosted_id)%colloid%conc_mob(:))* &
                1000.d0 ! convert kg water/L water -> kg water/m^3 water
          endif
        case(MASS_RATE_SS)
          Res(istartaq:iendaq) = -source_sink%tran_condition% &
                 cur_constraint_coupler%rt_auxvar%total(:,iphase) ! actually moles/sec
          if (reaction%ncoll > 0) then
            option%io_buffer = 'Need to implement MASS_RATE_SS source/sink ' // &
                               'term correctly for colloids'
            call printErrMsg(option)
            Res(istartcoll:iendcoll) = -source_sink%tran_condition% &
                 cur_constraint_coupler%rt_auxvar%colloid%conc_mob(:) ! actually moles/sec
          endif          
        case default
          if (qsrc > 0) then ! injection
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              Res(istartaq:iendaq) = -qsrc* & ! m^3 water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%total(:,iphase)*1000.d0
              if (reaction%ncoll > 0) then
                Res(istartcoll:iendcoll) = -qsrc* & ! m^3 water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%colloid%conc_mob(:)*1000.d0
              endif                     
            else ! mass
              Res(istartaq:iendaq) = -qsrc* & ! kg water / sec
                     source_sink%tran_condition%cur_constraint_coupler% &
                     rt_auxvar%total(:,iphase)/ &
                     global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                     1000.d0
              if (reaction%ncoll > 0) then  ! needs to be moles/sec
                Res(istartcoll:iendcoll) = -qsrc* & ! kg water / sec
                    source_sink%tran_condition%cur_constraint_coupler% &
                    rt_auxvar%colloid%conc_mob(:)/ &
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0
              endif                     
            endif
          else ! extraction
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              Res(istartaq:iendaq) = -qsrc*rt_aux_vars(ghosted_id)%total(:,iphase)*1000.d0
              if (reaction%ncoll > 0) then
                Res(istartcoll:iendcoll) = -qsrc* & ! m^3 water / sec
                    rt_aux_vars(ghosted_id)%colloid%conc_mob(:)*1000.d0
              endif               
            else
              Res(istartaq:iendaq) = -qsrc* &
                    rt_aux_vars(ghosted_id)%total(:,iphase)/ &
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
              if (reaction%ncoll > 0) then  ! needs to be moles/sec
                Res(istartcoll:iendcoll) = -qsrc* & ! kg water / sec
                    rt_aux_vars(ghosted_id)%colloid%conc_mob(:)/ & 
                    global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)* &
                    1000.d0 ! convert kg water/L water -> kg water/m^3 water
              endif                     
            endif
          endif
      end select
!      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
!        rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_aux_vars_ss(ghosted_id)%mass_balance_delta(:,iphase) + Res
!      endif
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo

#ifdef CHUAN_CO2
  select case(option%iflowmode)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      source_sink => patch%source_sinks%first 
      do 
        if (.not.associated(source_sink)) exit

        msrc(:) = source_sink%flow_condition%pressure%dataset%cur_value(:)
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
                if(abs(reaction%species_idx%co2_gas_id) == ieqgas) then
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
  use Level_module
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
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
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

  call PetscLogEventEnd(logging%event_rt_jacobian1,ierr)

  call PetscLogEventBegin(logging%event_rt_jacobian2,ierr)
  
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
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Jdn(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum

  PetscReal :: coef_up(realization%option%nphase)
  PetscReal :: coef_dn(realization%option%nphase)
    
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


  ! Get pointer to Vector data
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc, tor_loc_p, ierr)

  if (option%use_samr) then
#if 0
    ! again this should be broken
!!!!!!!!!!#ifndef REVISED_TRANSPORT
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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      call TFluxDerivativeAdv(rt_aux_vars(ghosted_id_up), &
                              global_aux_vars(ghosted_id_up), &
                              rt_aux_vars(ghosted_id_dn), &
                              global_aux_vars(ghosted_id_dn), &
                              cur_connection_set%area(iconn), &
                              rt_parameter,option, &
                              patch%internal_velocities(:,sum_connection),Jup,Jdn)

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

      if (patch%imat(ghosted_id) <= 0) cycle

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

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

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
                               patch%internal_velocities(:,sum_connection),Jup,Jdn)

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

      if (patch%imat(ghosted_id) <= 0) cycle

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
#endif  
!!!!!!!!!!!!#endif ! #ifndef REVISED_TRANSPORT
  else ! !use_samr

  ! Interior Flux Terms -----------------------------------
  ! must zero out Jacobian blocks

  call PetscLogEventBegin(logging%event_rt_jacobian_flux,ierr)

  Jup = 0.d0  
  Jdn = 0.d0  
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

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  call PetscLogEventEnd(logging%event_rt_jacobian_flux,ierr)
  
  ! Boundary Flux Terms -----------------------------------
  ! must zero out Jacobian block

  call PetscLogEventBegin(logging%event_rt_jacobian_fluxbc,ierr)

  Jdn = 0.d0
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
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  call PetscLogEventEnd(logging%event_rt_jacobian_fluxbc,ierr)  
  endif ! #else AMR_FLUX

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
  
  implicit none

  interface
     subroutine SAMRSetJacobianSrcCoeffsOnPatch(which_pc, p_application, p_patch) 
#include "finclude/petscsys.h"

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
  PetscInt :: iconn
  PetscReal :: qsrc, rdum
  PetscBool :: volumetric
  
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
  call PetscLogEventBegin(logging%event_rt_jacobian_accum,ierr)  
#if 1  
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RTAccumulationDerivative(rt_aux_vars(ghosted_id), &
                                    global_aux_vars(ghosted_id), &
                                    porosity_loc_p(ghosted_id), &
                                    volume_p(local_id),reaction,option,Jup) 
      if (reaction%neqsorb > 0 .and. reaction%kinmr_nrate <= 0) then
        call RAccumulationSorbDerivative(rt_aux_vars(ghosted_id), &
                                         global_aux_vars(ghosted_id), &
                                         volume_p(local_id),reaction, &
                                         option,Jup)
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

      if (patch%imat(ghosted_id) <= 0) cycle

      istartaq = reaction%offset_aq + 1
      iendaq = reaction%offset_aq + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_coll + 1
        iendcoll = reaction%offset_coll + reaction%ncoll
      endif
      
      Jup = 0.d0
      select case(source_sink%tran_condition%itype)
        case(EQUILIBRIUM_SS)
          do idof = istartaq, iendaq
            Jup(idof,idof) = 1.d-6* &
                                 porosity_loc_p(ghosted_id)* &
                                 global_aux_vars(ghosted_id)%sat(option%liquid_phase)* &
                                 volume_p(local_id)
          enddo
          if (reaction%ncoll > 0) then
             do idof = istartcoll, iendcoll
              Jup(idof,idof) = 1.d-6* &
                               porosity_loc_p(ghosted_id)* &
                               global_aux_vars(ghosted_id)%sat(option%liquid_phase)* &
                               volume_p(local_id)
            enddo
          endif
        case(MASS_RATE_SS)
        case default
          ! units = kg water/sec
          if (qsrc < 0) then ! extraction
            if (volumetric) then ! qsrc is volumetric; must be converted to mass
              ! qsrc = m^3 water/sec
              do idof = istartaq, iendaq
                ! m^3 water/sec * kg water/m^3 water = kg water/sec
                Jup(idof,idof) = -qsrc*global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)
              enddo
              if (reaction%ncoll > 0) then
                do idof = istartcoll, iendcoll
                  Jup(idof,idof) = -qsrc*global_aux_vars(ghosted_id)%den_kg(option%liquid_phase)
                enddo
              endif              
            else ! qsrc is mass -> kg water/sec
              do idof = istartaq, iendaq
                Jup(idof,idof) = -qsrc
              enddo
              if (reaction%ncoll > 0) then
                do idof = istartcoll, iendcoll
                  Jup(idof,idof) = -qsrc
                enddo
              endif                 
            endif
          endif
      end select
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr) 
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
      istartaq = offset + reaction%offset_aq + 1
      iendaq = offset + reaction%offset_aq + reaction%naqcomp
      if (reaction%ncoll > 0) then
        istartcoll = offset + reaction%offset_coll + 1
        iendcoll = offset + reaction%offset_coll + reaction%ncoll
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
  PetscBool :: update_bcs
  PetscBool :: update_activity_coefs
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      ! PETSC_TRUE to update cells
      call RTUpdateAuxVarsPatch(realization,PETSC_TRUE,update_bcs, &
                                update_activity_coefs)
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
  
  call GridVecGetArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)

  if (update_cells) then

    call PetscLogEventBegin(logging%event_rt_auxvars,ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials

      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*reaction%ncomp
      istartaq = offset + reaction%offset_aq + 1
      iendaq = offset + reaction%offset_aq + reaction%naqcomp
      
      patch%aux%RT%aux_vars(ghosted_id)%pri_molal = xx_loc_p(istartaq:iendaq)
      if (reaction%ncoll > 0) then
        istartcoll = offset + reaction%offset_coll + 1
        iendcoll = offset + reaction%offset_coll + reaction%ncoll
        patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_mob = xx_loc_p(istartcoll:iendcoll)* &
          patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)*1.d-3
      endif
      
      if (compute_activity_coefs) then
        call RActivityCoefficients(patch%aux%RT%aux_vars(ghosted_id), &
                                   patch%aux%Global%aux_vars(ghosted_id), &
                                   reaction,option)
        if(option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE)then
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

#ifdef DASVYAT
!		do iconn=1,6
! 			write(*,*) "total", iconn, patch%aux%RT%aux_vars_bc(iconn)%total(1,1)
!        end do 
#endif          

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
#ifdef DASVYAT
!		write(*,*) "basis_molarity_p",basis_molarity_p,"den_kg",patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(1)
#endif
        
        if (patch%imat(ghosted_id) <= 0) cycle

        offset = (ghosted_id-1)*reaction%ncomp
        istartaq_loc = reaction%offset_aq + 1
        iendaq_loc = reaction%offset_aq + reaction%naqcomp
        istartaq = offset + istartaq_loc
        iendaq = offset + iendaq_loc
    
        if (reaction%ncoll > 0) then
          istartcoll_loc = reaction%offset_coll + 1
          iendcoll_loc = reaction%offset_coll + reaction%ncoll
          istartcoll = offset + istartcoll_loc
          iendcoll = offset + iendcoll_loc
        endif

!       if (option%iflowmode /= MPH_MODE .or. icall>1) then
        if (option%iflowmode /= MPH_MODE .and. option%iflowmode /= FLASH2_MODE)then
!       Note: the  DIRICHLET_BC is not time dependent in this case (icall)    
        select case(boundary_condition%tran_condition%itype)
            case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
              ! since basis_molarity is in molarity, must convert to molality
                ! by dividing by density of water (mol/L -> mol/kg)
              xxbc(istartaq_loc:iendaq_loc) = basis_molarity_p(1:reaction%naqcomp) / &
                patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
              if (reaction%ncoll > 0) then
                xxbc(istartcoll_loc:iendcoll_loc) = basis_coll_conc_p(1:reaction%ncoll) / &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
              endif
            case(DIRICHLET_ZERO_GRADIENT_BC)
  !geh            do iphase = 1, option%nphase
                if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! same as dirichlet above
                  xxbc(istartaq_loc:iendaq_loc) = basis_molarity_p(1:reaction%naqcomp) / &
                    patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
                  if (reaction%ncoll > 0) then
                    xxbc(istartcoll_loc:iendcoll_loc) = basis_coll_conc_p(1:reaction%ncoll) / &
                      patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
                  endif
                else
                  ! same as zero_gradient below
                  xxbc(istartaq_loc:iendaq_loc) = xx_loc_p(istartaq:iendaq)
                  if (reaction%ncoll > 0) then
                    xxbc(istartcoll_loc:iendcoll_loc) = basis_coll_conc_p(1:reaction%ncoll) / &
                      patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
                  endif
                endif
  !geh          enddo
            case(ZERO_GRADIENT_BC)
              xxbc(istartaq_loc:iendaq_loc) = xx_loc_p(istartaq:iendaq)
              if (reaction%ncoll > 0) then
                xxbc(istartcoll_loc:iendcoll_loc) = basis_coll_conc_p(1:reaction%ncoll) / &
                  patch%aux%Global%aux_vars_bc(sum_connection)%den_kg(iphase) * 1000.d0
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
            if(option%iflowmode == MPH_MODE .or. option%iflowmode == FLASH2_MODE)then
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
            case(DIRICHLET_ZERO_GRADIENT_BC)
                if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! don't need to do anything as the constraint below provides all
                  ! the concentrations, etc.
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
              boundary_condition%tran_condition%cur_constraint_coupler%surface_complexes, &
              boundary_condition%tran_condition%cur_constraint_coupler%colloids, &
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
      enddo
      boundary_condition => boundary_condition%next
    enddo

    patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE

    call PetscLogEventEnd(logging%event_rt_auxvars_bc,ierr)

  endif 
  
  call GridVecRestoreArrayF90(grid,field%tran_xx_loc,xx_loc_p, ierr)
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


  if(.not.(option%use_samr)) then
     call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                        MPI_MAX,option%mycomm,ierr)
     
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
  character(len=2) :: free_mol_char, tot_mol_char
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  PetscInt :: i
  
  option => realization%option
  reaction => realization%reaction
  
  string = ''
  
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
  
  if (reaction%print_pH .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%h_ion_id > 0) then
      if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-pH"'')') icolumn
      else
        write(string2,'('',"pH"'')') 
      endif
      string = trim(string) // trim(string2)
    endif
  endif
  
  if (reaction%print_total_component) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_tot_'',a,''"'')') icolumn, &
            trim(reaction%primary_species_names(i)), trim(tot_mol_char)
        else
          write(string2,'('',"'',a,''_tot_'',a,''"'')') &
            trim(reaction%primary_species_names(i)), trim(tot_mol_char)
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  if (reaction%print_free_ion) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_free_'',a,''"'')') icolumn, &
            trim(reaction%primary_species_names(i)), trim(free_mol_char)
        else
          write(string2,'('',"'',a,''_free_'',a,''"'')') &
            trim(reaction%primary_species_names(i)), trim(free_mol_char)
        endif
        string = trim(string) // trim(string2)
      endif
    enddo  
  endif
    
  if (reaction%print_act_coefs) then
    do i=1,reaction%naqcomp
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
  
  do i=1,realization%reaction%neqsrfcplxrxn
    if (reaction%eqsrfcplx_site_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%eqsrfcplx_site_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%eqsrfcplx_site_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,realization%reaction%neqsrfcplx
    if (reaction%eqsrfcplx_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%eqsrfcplx_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%eqsrfcplx_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,realization%reaction%nkinsrfcplxrxn
    if (reaction%kinsrfcplx_site_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%kinsrfcplx_site_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%kinsrfcplx_site_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo
  
  do i=1,realization%reaction%nkinsrfcplx
    if (reaction%kinsrfcplx_print(i)) then
      if (icolumn > -1) then
        icolumn = icolumn + 1  
        write(string2,'('',"'',i2,''-'',a,''"'')') icolumn, &
          trim(reaction%kinsrfcplx_names(i))
      else
        write(string2,'('',"'',a,''"'')') trim(reaction%kinsrfcplx_names(i))
      endif
      string = trim(string) // trim(string2)
    endif
  enddo

  if (associated(reaction%kd_print)) then
    do i=1,reaction%naqcomp
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
    do i=1,reaction%naqcomp
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
  
  if (associated(reaction%total_sorb_mobile_print)) then
    do i=1,reaction%ncollcomp
      if (reaction%total_sorb_mobile_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_total_sorb_mob"'')') icolumn, &
            trim(reaction%colloid_species_names(i))
        else
          write(string2,'('',"'',a,''_total_sorb_mob"'')') &
            trim(reaction%colloid_species_names(i))
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  if (reaction%print_colloid) then
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_col_mob_'',a,''"'')') icolumn, &
            trim(reaction%colloid_names(i)), trim(tot_mol_char)
        else
          write(string2,'('',"'',a,''_col_mob_'',a,''"'')') &
            trim(reaction%colloid_names(i)), trim(tot_mol_char)
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-'',a,''_col_imb_'',a,''"'')') icolumn, &
            trim(reaction%colloid_names(i)), trim(tot_mol_char)
        else
          write(string2,'('',"'',a,''_col_imb_'',a,''"'')') &
            trim(reaction%colloid_names(i)), trim(tot_mol_char)
        endif
        string = trim(string) // trim(string2)
      endif
    enddo
  endif
  
  if (reaction%print_age) then
    if (reaction%species_idx%tracer_age_id > 0) then
        if (icolumn > -1) then
          icolumn = icolumn + 1
          write(string2,'('',"'',i2,''-Tracer_Age"'')') icolumn
        else
          write(string2,'('',"Tracer_Age"'')') 
        endif
        string = trim(string) // trim(string2)
    endif
  endif
    
  RTGetTecplotHeader = string

end function RTGetTecplotHeader

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
      call RTJumpStartKineticSorptionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
                              
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
  
  ! This subroutine assumes that the auxilliary variables are current!

  if (reaction%kinmr_nrate > 0) then
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
  use Level_module
  use Grid_module
  use Option_module
  use Field_module
  
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscBool :: checkpoint
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, pointer :: vec_p(:)

  PetscBool :: checkpoint_flag(realization%reaction%naqcomp)
  PetscInt :: i, j, irxn, icomp, icplx, ncomp, ncplx, irate
  PetscInt :: local_id
  PetscErrorCode :: ierr
  
  option => realization%option
  reaction => realization%reaction
  field => realization%field
  
  checkpoint_flag = PETSC_FALSE

  ! Loop over sorption reactions to find the necessary components
  
  do irxn = 1, reaction%neqsrfcplxrxn
    ncplx = reaction%eqsrfcplx_rxn_to_complex(0,irxn)
    do j = 1, ncplx
      icplx = reaction%eqsrfcplx_rxn_to_complex(j,irxn)
      ncomp = reaction%eqsrfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = reaction%eqsrfcplxspecid(i,icplx)
        checkpoint_flag(icomp) = PETSC_TRUE
      enddo
    enddo
  enddo

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      rt_auxvars => cur_patch%aux%RT%aux_vars
      grid => cur_patch%grid
      do icomp = 1, reaction%naqcomp
        if (checkpoint_flag(icomp)) then
          do irate = 1, reaction%kinmr_nrate
            if (checkpoint) then
              call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
              do local_id = 1, grid%nlmax
                vec_p(local_id) = &
                  rt_auxvars(grid%nL2G(local_id))%kinmr_total_sorb(icomp,irate)
              enddo
              call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
              call VecView(field%work,viewer,ierr)
            else
              call VecLoad(viewer,field%work,ierr)
              call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
              do local_id = 1, grid%nlmax
                rt_auxvars(grid%nL2G(local_id))%kinmr_total_sorb(icomp,irate) = &
                   vec_p(local_id)
              enddo
              call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
            endif
          enddo
        endif
      enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RTCheckpointKineticSorption

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
  use Option_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

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

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

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
      realization%patch => cur_patch
      call RTDestroyPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  
  ! taken care of in auxilliary.F90
  
end subroutine RTDestroyPatch

end module Reactive_Transport_module
