module Richards_module

  use Richards_Aux_module
  use Richards_Common_module
  use Global_Aux_module
  use Material_Aux_class
#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsResidual, &
         RichardsJacobian, &
         RichardsTimeCut,&
         RichardsSetup, &
         RichardsInitializeTimestep, &
         RichardsUpdateAuxVars, &
         RichardsMaxChange, &
         RichardsUpdateSolution, &
         RichardsComputeMassBalance, &
         RichardsDestroy, &
         RichardsCheckUpdatePre, &
         RichardsCheckUpdatePost, &
         RichardsUpdateSurfacePress

contains

! ************************************************************************** !

subroutine RichardsTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
 
  use Realization_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscViewer :: viewer

  option => realization%option
  field => realization%field

  
  call VecCopy(field%flow_yy,field%flow_xx,ierr)

  if (option%mimetic) then
    call VecCopy(field%flow_yy_faces, field%flow_xx_faces, ierr)
    call RichardsUpdateAuxVars(realization)
    call VecCopy(field%flow_yy,field%flow_xx,ierr)
!    read(*,*)    

  endif

  call RichardsInitializeTimestep(realization)  
 
end subroutine RichardsTimeCut

! ************************************************************************** !

subroutine RichardsSetup(realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  call RichardsSetupPatch(realization)
  call RichardsSetPlotVariables(realization)

end subroutine RichardsSetup

! ************************************************************************** !

subroutine RichardsSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Realization_class
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
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, ierr
  type(richards_auxvar_type), pointer :: rich_auxvars(:)  
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)  
  type(richards_auxvar_type), pointer :: rich_auxvars_ss(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Richards => RichardsAuxCreate()
  allocate(patch%aux%Richards%richards_parameter%sir(option%nphase, &
                                  size(patch%saturation_function_array)))
  do i = 1, size(patch%saturation_function_array)
    patch%aux%Richards%richards_parameter%sir(:,patch%saturation_function_array(i)%ptr%id) = &
      patch%saturation_function_array(i)%ptr%Sr(:)
  enddo
  
  ! allocate auxvar data structures for all grid cells  
  allocate(rich_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RichardsAuxVarInit(rich_auxvars(ghosted_id),option)
  enddo
  patch%aux%Richards%auxvars => rich_auxvars
  patch%aux%Richards%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  if (sum_connection > 0) then
    allocate(rich_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_auxvars_bc(iconn),option)
    enddo
    patch%aux%Richards%auxvars_bc => rich_auxvars_bc
  endif
  patch%aux%Richards%num_aux_bc = sum_connection
  
  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(rich_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_auxvars_ss(iconn),option)
    enddo
    patch%aux%Richards%auxvars_ss => rich_auxvars_ss
  endif
  patch%aux%Richards%num_aux_ss = sum_connection

#ifdef YE_FLUX
  allocate(patch%internal_fluxes(1,1,ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list)))
  patch%internal_fluxes = 0.d0
#endif

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RichardsCreateZeroArray(patch,option)


end subroutine RichardsSetupPatch

! ************************************************************************** !

subroutine RichardsCheckUpdatePre(line_search,P,dP,changed,realization,ierr)
  ! 
  ! Checks update prior to update
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/12
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  ! ignore changed flag for now.
  PetscBool :: changed
  type(realization_type) :: realization
  
  PetscReal, pointer :: P_p(:)
  PetscReal, pointer :: dP_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P_R, P0, P1, delP
  PetscReal :: scale, sat, sat_pert, pert, pc_pert, press_pert, delP_pert
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => realization%patch%aux%Richards%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars

  if (dabs(option%saturation_change_limit) > 0.d0) then

    patch => realization%patch

    call VecGetArrayF90(dP,dP_p,ierr)
    call VecGetArrayF90(P,P_p,ierr)

    pert =dabs(option%saturation_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      sat = global_auxvars(ghosted_id)%sat(1)
      sat_pert = sat - sign(1.d0,sat-0.5d0)*pert
      call SatFuncGetCapillaryPressure(pc_pert,sat_pert, &
             patch%saturation_function_array( &
               patch%sat_func_id(ghosted_id))%ptr,option)
      press_pert = option%reference_pressure - pc_pert
      P0 = P_p(local_id)
      delP = dP_p(local_id)
      delP_pert = dabs(P0 - press_pert)
      if (delP_pert < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),delP_pert,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),delP_pert),delP)
      dP_p(local_id) = delP
    enddo
    
    call VecRestoreArrayF90(dP,dP_p,ierr)
    call VecRestoreArrayF90(P,P_p,ierr)

  endif

  if (dabs(option%pressure_dampening_factor) > 0.d0) then
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = option%pressure_dampening_factor

    call VecGetArrayF90(dP,dP_p,ierr)
    call VecGetArrayF90(P,P_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
    do local_id = 1, grid%nlmax
      delP = dP_p(local_id)
      P0 = P_p(local_id)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1 
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dP_p(local_id) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dP,dP_p,ierr)
    call VecRestoreArrayF90(P,P_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
  endif

end subroutine RichardsCheckUpdatePre

! ************************************************************************** !

subroutine RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                                   P1_changed,realization,ierr)
  ! 
  ! Checks update after to update
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/12
  ! 

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  type(realization_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  
  PetscReal, pointer :: P1_p(:)
  PetscReal, pointer :: dP_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: Res(1)
  PetscReal :: inf_norm
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  rich_auxvars => realization%patch%aux%Richards%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  material_auxvars => realization%patch%aux%Material%auxvars
  
  dP_changed = PETSC_FALSE
  P1_changed = PETSC_FALSE
  
  if (option%check_stomp_norm) then
    call VecGetArrayF90(dP,dP_p,ierr)
    call VecGetArrayF90(P1,P1_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
    
    inf_norm = 0.d0
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
    
      call RichardsAccumulation(rich_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                option,Res)
      inf_norm = max(inf_norm,min(dabs(dP_p(local_id)/P1_p(local_id)), &
                                  dabs(r_p(local_id)/Res(1))))
    enddo
    call MPI_Allreduce(inf_norm,option%stomp_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    call VecRestoreArrayF90(dP,dP_p,ierr)
    call VecRestoreArrayF90(P1,P1_p,ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr)
  endif
  
end subroutine RichardsCheckUpdatePost

! ************************************************************************** !

subroutine RichardsComputeMassBalance(realization,mass_balance)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_class

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
  
  mass_balance = 0.d0
  
  call RichardsComputeMassBalancePatch(realization,mass_balance)

end subroutine RichardsComputeMassBalance

! ************************************************************************** !

subroutine RichardsComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_auxvars(ghosted_id)%den_kg* &
      global_auxvars(ghosted_id)%sat* &
      material_auxvars(ghosted_id)%porosity* &
      material_auxvars(ghosted_id)%volume
  enddo

end subroutine RichardsComputeMassBalancePatch

! ************************************************************************** !

subroutine RichardsZeroMassBalDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine RichardsZeroMassBalDeltaPatch

! ************************************************************************** !

subroutine RichardsUpdateMassBalancePatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance = &
      patch%aux%Global%auxvars(iconn)%mass_balance + &
      patch%aux%Global%auxvars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance = &
        global_auxvars_ss(iconn)%mass_balance + &
        global_auxvars_ss(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

end subroutine RichardsUpdateMassBalancePatch

! ************************************************************************** !

subroutine RichardsUpdatePermPatch(realization)
  ! 
  ! Updates the permeability based on pressure
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 01/09/12
  ! 

  use Grid_module
  use Realization_class
  use Option_module
  use Discretization_module
  use Patch_module
  use Field_module
  use Material_module
  use Material_Aux_class
  use Variables_module
  
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(discretization_type), pointer :: discretization
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscReal :: scale
  PetscReal :: p_min, p_max, permfactor_max
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  grid => patch%grid
  material_property_array => realization%material_property_array
  material_auxvars => patch%aux%Material%auxvars

  if (.not.associated(patch%imat)) then
    option%io_buffer = 'Materials IDs not present in run.  Material ' // &
      ' properties cannot be updated without material ids'
    call printErrMsg(option)
  endif
  
  call VecGetArrayF90(field%perm0_xx,perm0_xx_p,ierr)
  call VecGetArrayF90(field%perm0_zz,perm0_zz_p,ierr)
  call VecGetArrayF90(field%perm0_yy,perm0_yy_p,ierr)
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    p_min = material_property_array(patch%imat(ghosted_id))%ptr%min_pressure
    p_max = material_property_array(patch%imat(ghosted_id))%ptr%max_pressure
    permfactor_max = material_property_array(patch%imat(ghosted_id))%ptr% &
                     max_permfactor
    if (xx_loc_p(local_id) < p_min) then
      scale = 1
    else 
      if (xx_loc_p(local_id) < p_max) then
        scale = (xx_loc_p(local_id) - p_min)/ &
                (p_max - p_min)*(permfactor_max - 1.d0) + 1.d0
      else
        scale = permfactor_max
      endif
    endif
    !geh: this is a kludge for gfortran.  the code reports errors when 
    !     material_auxvars(ghosted_id)%permeability is used.
    !TODO(geh): test with Intel!
    perm_ptr => material_auxvars(ghosted_id)%permeability
    perm_ptr(perm_xx_index) = perm0_xx_p(local_id)*scale
    perm_ptr(perm_yy_index) = perm0_yy_p(local_id)*scale
    perm_ptr(perm_zz_index) = perm0_zz_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_xx_index) = &
!      perm0_xx_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_yy_index) = &
!      perm0_yy_p(local_id)*scale
!    material_auxvars(ghosted_id)%permeability(perm_zz_index) = &
!      perm0_zz_p(local_id)*scale
  enddo
  
  call VecRestoreArrayF90(field%perm0_xx,perm0_xx_p,ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm0_zz_p,ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm0_yy_p,ierr)
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)

  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

  
end subroutine RichardsUpdatePermPatch

! ************************************************************************** !

subroutine RichardsUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Richards problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Richards_MFD_module

  type(realization_type) :: realization
  
  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
    call RichardsUpdateAuxVarsPatchMFDLP(realization)
  else
    call RichardsUpdateAuxVarsPatch(realization)
  end if  

end subroutine RichardsUpdateAuxVars

! ************************************************************************** !

subroutine RichardsUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Richards problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  
  use Richards_LSM_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)
  type(richards_auxvar_type), pointer :: rich_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase, i
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  Vec :: phi
  
  call PetscLogEventBegin(logging%event_r_auxvars,ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  rich_auxvars_ss => patch%aux%Richards%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    call RichardsAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id),rich_auxvars(ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       option)   
  enddo

  call PetscLogEventEnd(logging%event_r_auxvars,ierr)

  call PetscLogEventBegin(logging%event_r_auxvars_bc,ierr)

  ! boundary conditions
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

      select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC,HET_SURF_SEEPAGE_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC)
          xxbc(1) = xx_loc_p(ghosted_id)
      end select
     
 
      call RichardsAuxVarCompute(xxbc(1),rich_auxvars_bc(sum_connection), &
                         global_auxvars_bc(sum_connection), &
                         material_auxvars(ghosted_id), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         option)
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

      call RichardsAuxVarCopy(rich_auxvars(ghosted_id), &
                              rich_auxvars_ss(sum_connection),option)
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_ss(sum_connection),option)

    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)

  ! Compute gradient using a least squares approach at each control volume
  if(realization%discretization%lsm_flux_method) then
    call RichardsUpdateLSMAuxVarsPatch(realization)
  endif

  patch%aux%Richards%auxvars_up_to_date = PETSC_TRUE

  call PetscLogEventEnd(logging%event_r_auxvars_bc,ierr)

end subroutine RichardsUpdateAuxVarsPatch

! ************************************************************************** !

subroutine RichardsInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 

  use Realization_class
  use Field_module 
  
  implicit none
  


  type(realization_type) :: realization

  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(field_type), pointer :: field


  field => realization%field


  call RichardsUpdateFixedAccum(realization)

!   call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy.out', &
!                              viewer,ierr)
!    call VecView(field%flow_xx_faces, viewer, ierr)
!    call VecView(field%flow_yy, viewer, ierr)
!
!    call PetscViewerDestroy(viewer,ierr)
!    write(*,*) "Flow_yy" 
!    read(*,*)    

end subroutine RichardsInitializeTimestep

! ************************************************************************** !

subroutine RichardsUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/08
  ! 

  use Realization_class
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization

  type(field_type), pointer :: field
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
 
  if (realization%option%mimetic) then 
     call VecCopy(field%flow_xx_faces, field%flow_yy_faces, ierr)  
  end if

  call VecCopy(field%flow_xx,field%flow_yy,ierr)   

  call RichardsUpdateSolutionPatch(realization)

end subroutine RichardsUpdateSolution

! ************************************************************************** !

subroutine RichardsUpdateSolutionPatch(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/13/08
  ! 

  use Realization_class
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call RichardsUpdateMassBalancePatch(realization)
  endif
  
  if (realization%option%update_flow_perm) then
    call RichardsUpdatePermPatch(realization)
  endif

end subroutine RichardsUpdateSolutionPatch

! ************************************************************************** !

subroutine RichardsUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class

  type(realization_type) :: realization
  
  call RichardsUpdateFixedAccumPatch(realization)

end subroutine RichardsUpdateFixedAccum

! ************************************************************************** !

subroutine RichardsUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, numfaces, jface, ghost_face_id, j
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx,xx_p, ierr)

!  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

!  numfaces = 6     ! hex only
!  allocate(sq_faces(numfaces))
!  allocate(faces_pr(numfaces))

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsAuxVarCompute(xx_p(local_id:local_id), &
                   rich_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                   material_auxvars(ghosted_id), &
                   patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                   option)
    call RichardsAccumulation(rich_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              option,accum_p(local_id:local_id))
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)


!  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)

#if 0
!  call RichardsNumericalJacTest(field%flow_xx,realization)
#endif

end subroutine RichardsUpdateFixedAccumPatch

! ************************************************************************** !

subroutine RichardsNumericalJacTest(xx,realization)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Realization_class
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
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  
  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr)
  call VecDuplicate(xx,res,ierr)
  call VecDuplicate(xx,res_pert,ierr)
  
  call MatCreate(option%mycomm,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof,grid%nlmax*option%nflowdof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call RichardsResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
     idof = icell
!    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call veccopy(xx,xx_pert,ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr)
      call RichardsResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr)
!    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine RichardsNumericalJacTest

! ************************************************************************** !

subroutine RichardsResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Field_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Mass_Transfer_module, only : mass_transfer_type
  use Material_module
  use Material_Aux_class
  use Variables_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(mass_transfer_type), pointer :: cur_mass_transfer
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

  ! pass #1 for internal and boundary flux terms
  call RichardsResidualPatch1(snes,xx,r,realization,ierr)

  ! pass #2 for everything else
  call RichardsResidualPatch2(snes,xx,r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  call PetscLogEventEnd(logging%event_r_residual,ierr)

  ! Mass Transfer
  if (associated(realization%flow_mass_transfer_list)) then
    cur_mass_transfer => realization%flow_mass_transfer_list
    do
      if (.not.associated(cur_mass_transfer)) exit
      call VecStrideScatter(cur_mass_transfer%vec,cur_mass_transfer%idof-1, &
                            r,ADD_VALUES,ierr)
      cur_mass_transfer => cur_mass_transfer%next
    enddo
  endif

end subroutine RichardsResidual

! ************************************************************************** !

subroutine RichardsResidualPatch1(snes,xx,r,realization,ierr)
  ! 
  ! Computes the interior flux and boundary flux
  ! terms of the residual equation on a single patch
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  use Richards_LSM_module
  
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
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof), v_darcy


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscViewer :: viewer
  PetscInt, pointer :: cell_neighbors(:,:)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call RichardsUpdateAuxVarsPatch(realization)
  patch%aux%Richards%auxvars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  patch%aux%Richards%auxvars_cell_pressures_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call RichardsZeroMassBalDeltaPatch(realization)
  endif

  select case(grid%itype)
    case(STRUCTURED_GRID)
      cell_neighbors => grid%structured_grid%cell_neighbors
    case(UNSTRUCTURED_GRID)
      option%io_buffer='GridComputeMinv() not implemented for unstructured grid.'
      call printErrMsg(option)
  end select

!  write(*,*) "RichardsResidual"
!  read(*,*)
! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr)

  r_p = 0.d0

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

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      select case (realization%discretization%hydr_flux_method)
        case (TWO_POINT_FLUX)
          call RichardsFlux(rich_auxvars(ghosted_id_up), &
                            global_auxvars(ghosted_id_up), &
                            material_auxvars(ghosted_id_up), &
                            richards_parameter%sir(1,icap_up), &
                            rich_auxvars(ghosted_id_dn), &
                            global_auxvars(ghosted_id_dn), &
                            material_auxvars(ghosted_id_dn), &
                            richards_parameter%sir(1,icap_dn), &
                            cur_connection_set%area(iconn), &
                            cur_connection_set%dist(:,iconn), &
                            option,v_darcy,Res)
        case (LSM_FLUX)
#if 0         
          option%io_buffer = 'RicardsLSM needs to be implemented with ' // &
                             'new material_aux_type.'
          call printErrMsg(option)
          call RichardsLSMFlux(rich_auxvars(ghosted_id_up), &
                               global_auxvars(ghosted_id_up), &
                               porosity_loc_p(ghosted_id_up), &
                               richards_parameter%sir(1,icap_up), &
                               dd_up,perm_up, &
                               rich_auxvars(ghosted_id_dn), &
                               global_auxvars(ghosted_id_dn), &
                               porosity_loc_p(ghosted_id_dn), &
                               richards_parameter%sir(1,icap_dn), &
                               dd_dn,perm_dn, &
                               cur_connection_set%area(iconn), &
                               cur_connection_set%dist(1:3,iconn), &
                               distance_gravity, &
                               upweight,distance,option,&
                               ghosted_id_up, ghosted_id_dn, &
                               cell_neighbors, &
                               grid%bnd_cell, &
                               v_darcy,Res)
#endif          
        case default
          option%io_buffer = 'Unknown hydr_flux_method '
          call printErrMsg(option)
      end select

      patch%internal_velocities(1,sum_connection) = v_darcy

#ifdef COMPUTE_INTERNAL_MASS_FLUX
      global_auxvars(local_id_up)%mass_balance_delta(1,1) = &
        global_auxvars(local_id_up)%mass_balance_delta(1,1) - Res(1)
#endif

#ifdef YE_FLUX
      patch%internal_fluxes(RICHARDS_PRESSURE_DOF,1,sum_connection) = Res(1)
#endif
#ifdef STORE_FLOWRATES
      patch%internal_fluxes(RICHARDS_PRESSURE_DOF,1,sum_connection) = Res(1)
#endif
      if (local_id_up>0) then
#ifdef PM_RICHARDS_DEBUG
  print *, 'Res interior up', local_id_up
  print *, Res(1)
#endif  
        r_p(local_id_up) = r_p(local_id_up) + Res(1)
      endif
         
      if (local_id_dn>0) then
#ifdef PM_RICHARDS_DEBUG
  print *, 'Res interior dn', local_id_dn
  print *, Res(1)
#endif  
        r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call RichardsBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                rich_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                rich_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                richards_parameter%sir(1,icap_dn), &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(:,iconn), &
                                option, &
                                v_darcy,Res)
      patch%boundary_velocities(1,sum_connection) = v_darcy
#ifdef STORE_FLOWRATES
      patch%boundary_fluxes(1,1,sum_connection) = Res(1)
#endif

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
        ! contribution to internal 
!        global_auxvars(ghosted_id)%mass_balance_delta(1) = &
!          global_auxvars(ghosted_id)%mass_balance_delta(1) + Res(1)
      endif

#ifdef PM_RICHARDS_DEBUG
  print *, 'Res bc', local_id
  print *, Res(1)
#endif

      r_p(local_id)= r_p(local_id) - Res(1)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(r, r_p, ierr)

  !read(*,*) local_id

end subroutine RichardsResidualPatch1

! ************************************************************************** !

subroutine RichardsResidualPatch2(snes,xx,r,realization,ierr)
  ! 
  ! Computes the accumulation and source/sink terms of
  ! the residual equation on a single patch
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  

  use Connection_module
  use Realization_class
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
  PetscViewer :: viewer
  PetscInt :: i
  PetscInt :: local_id, ghosted_id

  PetscReal, pointer :: r_p(:), accum_p(:)
  PetscReal :: qsrc, qsrc_mol
  PetscReal :: Res(realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal, pointer :: mmsrc(:)
  PetscReal, allocatable :: msrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: well_inj_water
  PetscReal :: Dq, dphi, v_darcy, ukvr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_ss => patch%aux%Richards%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%steady_state) then
    r_p = r_p - accum_p

    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RichardsAccumulation(rich_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                option,Res) 
#ifdef PM_RICHARDS_DEBUG
  print *, 'Res accum', local_id
  print *, Res(1)      
#endif
      r_p(local_id) = r_p(local_id) + Res(1)
    enddo
  endif

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

      if (source_sink%flow_condition%itype(1)/=HET_VOL_RATE_SS .and. &
          source_sink%flow_condition%itype(1)/=HET_MASS_RATE_SS .and. &
          source_sink%flow_condition%itype(1)/=WELL_SS) &
        qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O ! kg/sec -> kmol/sec
        case(SCALED_MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O* & ! kg/sec -> kmol/sec
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_VOL_RATE_SS)
          ! qsrc1 = m^3/sec
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* & ! flow = m^3/s
                     global_auxvars(ghosted_id)%den(1)                  ! den  = kmol/m^3
        case(HET_MASS_RATE_SS)
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O ! kg/sec -> kmol/sec
      
        case(WELL_SS) ! production well, SK 12/19/13
        ! if node pessure is lower than the given extraction pressure, shut it down
!  well parameter explanation
!   1. well status. 1 injection; -1 production; 0 shut in!
!   2. well factor [m^3],  the effective permeability [m^2/s]
!   3. bottomhole pressure:  [Pa]
!   4. max pressure: [Pa]
!   5. min pressure: [Pa]   
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)
    
        ! production well (well status = -1)
          if(dabs(well_status + 1.D0) < 1.D-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor 
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.D0) then ! outflow only
                ukvr = rich_auxvars(ghosted_id)%kvr
                if (ukvr < 1.e-20) ukvr = 0.D0
                v_darcy = 0.D0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  qsrc_mol = -1.d0*v_darcy*global_auxvars(ghosted_id)%den(1)
                endif
              endif
            endif
          endif 


      end select
          
      
      if (option%compute_mass_balance_new) then
        ! need to added global auxvar for src/sink
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) - &
          qsrc_mol
      endif
      r_p(local_id) = r_p(local_id) - qsrc_mol
      ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
      patch%ss_fluid_fluxes(1,sum_connection) = qsrc_mol / &
                                             global_auxvars(ghosted_id)%den(1)
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%Richards%inactive_cells_exist) then
    do i=1,patch%aux%Richards%n_zero_rows
      r_p(patch%aux%Richards%zero_rows_local(i)) = 0.d0
    enddo
  endif

#ifdef PM_RICHARDS_DEBUG
  print *, 'Residual'
  print *, r_p(:)
#endif

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  
end subroutine RichardsResidualPatch2

! ************************************************************************** !

subroutine RichardsJacobian(snes,xx,A,B,flag,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
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
  type(option_type), pointer :: option
  PetscReal :: norm
  
  call PetscLogEventBegin(logging%event_r_jacobian,ierr)

  option => realization%option

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
  call RichardsJacobianPatch1(snes,xx,J,J,flag,realization,ierr)

  ! pass #2 for everything else
  call RichardsJacobianPatch2(snes,xx,J,J,flag,realization,ierr)

  if (realization%debug%matview_Jacobian) then
#if 1  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr)
#else
!    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
!                               FILE_MODE_WRITE,viewer,ierr)
#endif    
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

#if 0
    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_dxx,viewer,ierr) 

    call PetscViewerDestroy(viewer,ierr)
 

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_yy,viewer,ierr) 
    call PetscViewerDestroy(viewer,ierr)
#endif

  call PetscLogEventEnd(logging%event_r_jacobian,ierr)
!  call printErrMsg(option)

  
end subroutine RichardsJacobian

! ************************************************************************** !

subroutine RichardsJacobianPatch1(snes,xx,A,B,flag,realization,ierr)
  ! 
  ! Computes the interior flux and boundary flux
  ! terms of the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
       
  

  use Connection_module
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  use Richards_LSM_module  
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscInt :: icap_up,icap_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt, pointer :: cell_neighbors(:,:)
  
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (associated(patch%aux%Richards%matrix_buffer)) then
      call MatrixBufferZero(patch%aux%Richards%matrix_buffer)
    else
      patch%aux%Richards%matrix_buffer => MatrixBufferCreate()
      call MatrixBufferInit(A,patch%aux%Richards%matrix_buffer,grid)
    endif
  endif
#endif

  select case(grid%itype)
    case(STRUCTURED_GRID)
      cell_neighbors => grid%structured_grid%cell_neighbors
    case(UNSTRUCTURED_GRID)
      option%io_buffer='GridComputeMinv() not implemented for unstructured grid.'
      call printErrMsg(option)
  end select

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

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      select case (realization%discretization%hydr_flux_method)
        case (TWO_POINT_FLUX)
          call RichardsFluxDerivative(rich_auxvars(ghosted_id_up), &
                                      global_auxvars(ghosted_id_up), &
                                      material_auxvars(ghosted_id_up), &
                                      richards_parameter%sir(1,icap_up), &
                                      rich_auxvars(ghosted_id_dn), &
                                      global_auxvars(ghosted_id_dn), &
                                      material_auxvars(ghosted_id_dn), &
                                      richards_parameter%sir(1,icap_dn), &
                                      cur_connection_set%area(iconn), &
                                      cur_connection_set%dist(-1:3,iconn),&
                                      option,&
                                      patch%saturation_function_array(icap_up)%ptr,&
                                      patch%saturation_function_array(icap_dn)%ptr,&
                                      Jup,Jdn)
        case (LSM_FLUX)
#if 0        
          call RichardsLSMFluxDerivative(rich_auxvars(ghosted_id_up), &
                                         global_auxvars(ghosted_id_up), &
                                         porosity_loc_p(ghosted_id_up), &
                                         richards_parameter%sir(1,icap_up), &
                                         dd_up,perm_up, &
                                         rich_auxvars(ghosted_id_dn), &
                                         global_auxvars(ghosted_id_dn), &
                                         porosity_loc_p(ghosted_id_dn), &
                                         richards_parameter%sir(1,icap_dn), &
                                         dd_dn,perm_dn, &
                                         cur_connection_set%area(iconn), &
                                         cur_connection_set%dist(1:3,iconn), &
                                         distance_gravity, &
                                         upweight, &
                                         distance, &
                                         option,&
                                         patch%saturation_function_array(icap_up)%ptr,&
                                         patch%saturation_function_array(icap_dn)%ptr,&
                                         grid%jacfac, ghosted_id_up, ghosted_id_dn, &
                                         cell_neighbors, &
                                         grid%x,grid%y,grid%z, &
                                         grid%bnd_cell, &
                                         Jup,Jdn)
#endif          
        case default
          option%io_buffer = 'Unknown hydr_flux_method '
          call printErrMsg(option)
      end select

      if (local_id_up > 0) then
#ifdef PM_RICHARDS_DEBUG
  print *, 'Jac up', local_id_up
  print *, Jup(1,1)
#endif      
#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id_up, &
                               ghosted_id_up,Jup(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id_up, &
                               ghosted_id_dn,Jdn(1,1))
        else
#endif
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                        Jup,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                        Jdn,ADD_VALUES,ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif
      if (local_id_dn > 0) then
#ifdef PM_RICHARDS_DEBUG
  print *, 'Jac dn', local_id_dn
  print *, Jdn(1,1)
#endif        
        Jup = -Jup
        Jdn = -Jdn
#ifdef BUFFER_MATRIX
        if (option%use_matrix_buffer) then
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id_dn, &
                               ghosted_id_dn,Jdn(1,1))
          call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id_dn, &
                               ghosted_id_up,Jup(1,1))
        else
#endif
          call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                        Jdn,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                        Jup,ADD_VALUES,ierr)
#ifdef BUFFER_MATRIX
        endif
#endif
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(option%mycomm,'jacobian_flux.bin',FILE_MODE_WRITE,viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id) 

      call RichardsBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                rich_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                rich_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                richards_parameter%sir(1,icap_dn), &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(:,iconn), &
                                option, &
                                patch%saturation_function_array(icap_dn)%ptr,&
                                Jdn)
      Jdn = -Jdn
#ifdef PM_RICHARDS_DEBUG
  print *, 'Jac dn bc', local_id
  print *, Jdn(1,1)      
#endif
#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jdn(1,1))
      else
#endif
        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                               ADD_VALUES,ierr)
#ifdef BUFFER_MATRIX
      endif
#endif
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(option%mycomm,'jacobian_bcflux.bin',FILE_MODE_WRITE,viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
end subroutine RichardsJacobianPatch1

! ************************************************************************** !

subroutine RichardsJacobianPatch2(snes,xx,A,B,flag,realization,ierr)
  ! 
  ! Computes the accumulation and source/sink terms of
  ! the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 
       
  

  use Connection_module
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal :: qsrc
  PetscInt :: icap
  PetscInt :: local_id, ghosted_id
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: flow_pc
  PetscViewer :: viewer
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: ukvr, Dq, dphi, v_darcy
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  if (.not.option%steady_state) then
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    icap = patch%sat_func_id(ghosted_id)
    call RichardsAccumDerivative(rich_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              option, &
                              patch%saturation_function_array(icap)%ptr,&
                              Jup) 
#ifdef PM_RICHARDS_DEBUG
  print *, 'Jac accum'
  print *, Jup(1,1)     
#endif
#ifdef BUFFER_MATRIX
    if (option%use_matrix_buffer) then
      call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                           ghosted_id,Jup(1,1))
    else
#endif
      call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                             ADD_VALUES,ierr)
#ifdef BUFFER_MATRIX
    endif
#endif
  enddo
#endif
  endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_accum.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(option%mycomm,'jacobian_accum.bin',FILE_MODE_WRITE,viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    if(source_sink%flow_condition%itype(1)/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%itype(1)/=HET_MASS_RATE_SS .and. &
       source_sink%flow_condition%itype(1)/=WELL_SS) &
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%itype(1))
        case(MASS_RATE_SS,SCALED_MASS_RATE_SS,HET_MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_auxvars(ghosted_id)%dden_dp*FMWH2O* &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_VOL_RATE_SS)
          Jup(1,1) = -source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* &
                    rich_auxvars(ghosted_id)%dden_dp*FMWH2O
        case(WELL_SS) ! production well, SK 12/19/13
        ! if node pessure is lower than the given extraction pressure, shut it down
!  well parameter explanation
!   1. well status. 1 injection; -1 production; 0 shut in!
!   2. well factor [m^3],  the effective permeability [m^2/s]
!   3. bottomhole pressure:  [Pa]
!   4. max pressure: [Pa]
!   5. min pressure: [Pa]   
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)
    
        ! production well (well status = -1)
          if(dabs(well_status + 1.D0) < 1.D-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor 
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.D0) then ! outflow only
                ukvr = rich_auxvars(ghosted_id)%kvr
                if (ukvr < 1.e-20) ukvr = 0.D0
                v_darcy = 0.D0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  Jup(1,1) = -1.d0*(-Dq*rich_auxvars(ghosted_id)%dkvr_dp*dphi* &
                             global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*1.d0*global_auxvars(ghosted_id)%den(1) &
                             -Dq*ukvr*dphi*rich_auxvars(ghosted_id)%dden_dp)
                endif
              endif
            endif
          endif 
      end select
#ifdef BUFFER_MATRIX
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jup(1,1))
      else
#endif
        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)  
#ifdef BUFFER_MATRIX
      endif
#endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(option%mycomm,'jacobian_srcsink.bin',FILE_MODE_WRITE,viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
#ifdef BUFFER_MATRIX
  if (option%use_matrix_buffer) then
    if (patch%aux%Richards%inactive_cells_exist) then
      call MatrixBufferZeroRows(patch%aux%Richards%matrix_buffer, &
                                patch%aux%Richards%n_zero_rows, &
                                patch%aux%Richards%zero_rows_local_ghosted)
    endif
    call MatrixBufferSetValues(A,patch%aux%Richards%matrix_buffer)
  endif
#endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef BUFFER_MATRIX
  if (.not.option%use_matrix_buffer) then
#endif
    if (patch%aux%Richards%inactive_cells_exist) then
      qsrc = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%Richards%n_zero_rows, &
                            patch%aux%Richards%zero_rows_local_ghosted, &
                            qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
    endif
#ifdef BUFFER_MATRIX
  endif
#endif

end subroutine RichardsJacobianPatch2

! ************************************************************************** !

subroutine RichardsCreateZeroArray(patch,option)
  ! 
  ! Computes the zeroed rows for inactive grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  PetscInt :: flag
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr
    
  flag = 0
  grid => patch%grid
  
  n_zero_rows = 0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      n_zero_rows = n_zero_rows + option%nflowdof
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
      do idof = 1, option%nflowdof
        ncount = ncount + 1
        zero_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
        zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof+idof-1
      enddo
    endif
  enddo

  patch%aux%Richards%zero_rows_local => zero_rows_local
  patch%aux%Richards%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%Richards%n_zero_rows = n_zero_rows
  
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_MAX,option%mycomm,ierr)
  if (flag > 0) patch%aux%Richards%inactive_cells_exist = PETSC_TRUE
     
  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine RichardsCreateZeroArray

! ************************************************************************** !

subroutine RichardsMaxChange(realization)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/08
  ! 

  use Realization_Base_class
  use Option_module
  use Field_module
  
  implicit none
  
  class(realization_base_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  option => realization%option
  field => realization%field

  if (option%mimetic) then

    call VecWAXPY(field%flow_dxx_faces,-1.d0,field%flow_xx_faces,field%flow_yy_faces,ierr)
    call VecStrideNorm(field%flow_dxx_faces,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

    call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
    call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

  else

     call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
     call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

  end if

end subroutine RichardsMaxChange

! ************************************************************************** !

subroutine RichardsSetPlotVariables(realization)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Realization_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => realization%output_option%output_variable_list

  if (associated(list%first)) then
    return
  endif
  
  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
end subroutine RichardsSetPlotVariables

! ************************************************************************** !

subroutine RichardsPrintAuxVars(richards_auxvar,global_auxvar,cell_id)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/21/12
  ! 

  use Global_Aux_module

  implicit none

  type(richards_auxvar_type) :: richards_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: cell_id

  print *, '      cell: ', cell_id
  print *, '  pressure: ', global_auxvar%pres(1)
  print *, 'saturation: ', global_auxvar%sat(1)
  print *, '   density: ', global_auxvar%den_kg(1)
  print *, '        pc: ', richards_auxvar%pc
  print *, '       kvr: ', richards_auxvar%kvr
  print *, '   dkvr_dp: ', richards_auxvar%dkvr_dp
  print *, '   dsat_dp: ', richards_auxvar%dsat_dp
  print *, '   dden_dp: ', richards_auxvar%dden_dp

end subroutine RichardsPrintAuxVars

! ************************************************************************** !

subroutine RichardsUpdateSurfacePress(realization)
  ! 
  ! This routine updates the boundary pressure condition corresponding on
  ! the top surface of the subsurface domain accounting for the amount of
  ! infilitration/exfiltration in the previous subsurface timestep.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/31/13
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  use String_module
  use EOS_Water_module

  implicit none

  type(realization_type) :: realization

#ifdef SURFACE_FLOW
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscReal :: den
  PetscReal :: surfpress_old
  PetscReal :: surfpress_new
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
    

  call EOSWaterdensity(option%reference_temperature,option%reference_pressure,den)

  ! boundary conditions
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if(StringCompare(boundary_condition%name,'from_surface_bc')) then

      if(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF) /= &
         HET_SURF_SEEPAGE_BC) then
        call printErrMsg(option,'from_surface_bc is not of type ' // &
                        'HET_SURF_SEEPAGE_BC')
      endif

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        surfpress_old = &
          boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)

        surfpress_new = surfpress_old - &
          patch%boundary_velocities(1,sum_connection)*option%flow_dt* &
          (abs(option%gravity(3)))*den

        surfpress_new = max(surfpress_new,option%reference_pressure)

        boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn) = &
          surfpress_new
      enddo
        
    else
      sum_connection = sum_connection + cur_connection_set%num_connections
    endif
 
    boundary_condition => boundary_condition%next

  enddo
#endif

end subroutine RichardsUpdateSurfacePress

! ************************************************************************** !

subroutine RichardsDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_class
  
  implicit none

  type(realization_type) :: realization
  
  call RichardsDestroyPatch(realization)

end subroutine RichardsDestroy

! ************************************************************************** !

subroutine RichardsDestroyPatch(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/09
  ! 

  use Realization_class

  implicit none

  type(realization_type) :: realization
  
  ! taken care of in auxiliary.F90

end subroutine RichardsDestroyPatch

end module Richards_module
