module Richards_module

  use Richards_Aux_module
  use Global_Aux_module
#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif
  
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


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsResidual,RichardsJacobian, &
         RichardsUpdateFixedAccum,RichardsTimeCut,&
         RichardsSetup, RichardsNumericalJacTest, &
         RichardsInitializeTimestep, RichardsUpdateAuxVars, &
         RichardsMaxChange, RichardsUpdateSolution, &
         RichardsGetTecplotHeader, RichardsComputeMassBalance, &
         RichardsDestroy, RichardsResidualMFD, RichardsJacobianMFD, &
!          RichardsInitialPressureReconstruction, &
         RichardsResidualMFDLP, RichardsJacobianMFDLP, &
         RichardsCheckUpdatePre, RichardsCheckUpdatePost

contains

! ************************************************************************** !
!
! RichardsTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsTimeCut(realization)
 
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
!
! RichardsSetup: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsSetup(realization)

  use Realization_class
  use Patch_module

  type(realization_type) :: realization
  
  call RichardsSetupPatch(realization)
  call RichardsSetPlotVariables(realization)

end subroutine RichardsSetup

! ************************************************************************** !
!
! RichardsSetupPatch: Creates arrays for auxiliary variables
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsSetupPatch(realization)

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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)  
  type(richards_auxvar_type), pointer :: rich_aux_vars_bc(:)  
  type(richards_auxvar_type), pointer :: rich_aux_vars_ss(:)  
  
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
  
  ! allocate aux_var data structures for all grid cells  
  allocate(rich_aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RichardsAuxVarInit(rich_aux_vars(ghosted_id),option)
  enddo
  patch%aux%Richards%aux_vars => rich_aux_vars
  patch%aux%Richards%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! aux_var data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  if (sum_connection > 0) then
    allocate(rich_aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_aux_vars_bc(iconn),option)
    enddo
    patch%aux%Richards%aux_vars_bc => rich_aux_vars_bc
  endif
  patch%aux%Richards%num_aux_bc = sum_connection
  
  ! count the number of source/sink connections and allocate
  ! aux_var data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(rich_aux_vars_ss(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_aux_vars_ss(iconn),option)
    enddo
    patch%aux%Richards%aux_vars_ss => rich_aux_vars_ss
  endif
  patch%aux%Richards%num_aux_ss = sum_connection
  
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RichardsCreateZeroArray(patch,option)


end subroutine RichardsSetupPatch

! ************************************************************************** !
!
! RichardsCheckUpdatePre: Checks update prior to update
! author: Glenn Hammond
! date: 02/13/12
!
! ************************************************************************** !
#ifndef HAVE_SNES_API_3_2
subroutine RichardsCheckUpdatePre(line_search,P,dP,changed,realization,ierr)
#else
subroutine RichardsCheckUpdatePre(snes_,P,dP,realization,changed,ierr)
#endif

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
 
  implicit none
  
#ifndef HAVE_SNES_API_3_2
  SNESLineSearch :: line_search
#else
  SNES :: snes_
#endif
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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P_R, P0, P1, delP
  PetscReal :: scale, sat, sat_pert, pert, pc_pert, press_pert, delP_pert
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  rich_aux_vars => realization%patch%aux%Richards%aux_vars
  global_aux_vars => realization%patch%aux%Global%aux_vars

  if (dabs(option%saturation_change_limit) > 0.d0) then

    patch => realization%patch

    call GridVecGetArrayF90(grid,dP,dP_p,ierr)
    call GridVecGetArrayF90(grid,P,P_p,ierr)

    pert =dabs(option%saturation_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      sat = global_aux_vars(ghosted_id)%sat(1)
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
    
    call GridVecRestoreArrayF90(grid,dP,dP_p,ierr)
    call GridVecRestoreArrayF90(grid,P,P_p,ierr)

  endif

  if (dabs(option%pressure_dampening_factor) > 0.d0) then
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = option%pressure_dampening_factor

    call GridVecGetArrayF90(grid,dP,dP_p,ierr)
    call GridVecGetArrayF90(grid,P,P_p,ierr)
    call GridVecGetArrayF90(grid,field%flow_r,r_p,ierr)
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
        call RichardsPrintAuxVars(rich_aux_vars(ghosted_id), &
                                  global_aux_vars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_aux_vars(ghosted_id), &
                                  global_aux_vars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dP_p(local_id) = scale*delP
      endif
    enddo
    call GridVecRestoreArrayF90(grid,dP,dP_p,ierr)
    call GridVecRestoreArrayF90(grid,P,P_p,ierr)
    call GridVecGetArrayF90(grid,field%flow_r,r_p,ierr)
  endif

end subroutine RichardsCheckUpdatePre

! ************************************************************************** !
!
! RichardsCheckUpdatePost: Checks update after to update
! author: Glenn Hammond
! date: 02/13/12
!
! ************************************************************************** !
#ifndef HAVE_SNES_API_3_2
subroutine RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                                   P1_changed,realization,ierr)
#else
subroutine RichardsCheckUpdatePost(snes_,P0,dP,P1,realization,dP_changed, &
                                   P1_changed,ierr)
#endif

  use Realization_class
  use Grid_module
  use Field_module
  use Option_module
 
  implicit none
  
#ifndef HAVE_SNES_API_3_2
  SNESLineSearch :: line_search
#else
  SNES :: snes_
#endif
  Vec :: P0
  Vec :: dP
  Vec :: P1
  type(realization_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  
  PetscReal, pointer :: P1_p(:)
  PetscReal, pointer :: dP_p(:)
  PetscReal, pointer :: volume_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: Res(1)
  PetscReal :: inf_norm
  PetscErrorCode :: ierr
  
  grid => realization%patch%grid
  option => realization%option
  field => realization%field
  rich_aux_vars => realization%patch%aux%Richards%aux_vars
  global_aux_vars => realization%patch%aux%Global%aux_vars
  
  dP_changed = PETSC_FALSE
  P1_changed = PETSC_FALSE
  
  if (option%check_stomp_norm) then
    call GridVecGetArrayF90(grid,dP,dP_p,ierr)
    call GridVecGetArrayF90(grid,P1,P1_p,ierr)
    call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
    call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%flow_r,r_p,ierr)
    
    inf_norm = 0.d0
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
    
      call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option,Res)
      inf_norm = max(inf_norm,min(dabs(dP_p(local_id)/P1_p(local_id)), &
                                  dabs(r_p(local_id)/Res(1))))
    enddo
    call MPI_Allreduce(inf_norm,option%stomp_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    call GridVecRestoreArrayF90(grid,dP,dP_p,ierr)
    call GridVecRestoreArrayF90(grid,P1,P1_p,ierr)
    call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
    call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%flow_r,r_p,ierr)
  endif
  
end subroutine RichardsCheckUpdatePost

! ************************************************************************** !
!
! RichardsComputeMassBalance: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsComputeMassBalance(realization,mass_balance)

  use Realization_class

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
  
  mass_balance = 0.d0
  
  call RichardsComputeMassBalancePatch(realization,mass_balance)

end subroutine RichardsComputeMassBalance

! ************************************************************************** !
!
! RichardsComputeMassBalancePatch: Initializes mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RichardsComputeMassBalancePatch(realization,mass_balance)
 
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
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  PetscReal, pointer :: volume_p(:), porosity_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_aux_vars => patch%aux%Global%aux_vars

  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density
    mass_balance = mass_balance + &
      global_aux_vars(ghosted_id)%den_kg* &
      global_aux_vars(ghosted_id)%sat* &
      porosity_loc_p(ghosted_id)*volume_p(local_id)
  enddo

  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine RichardsComputeMassBalancePatch


subroutine RichardsCheckMassBalancePatch(realization)

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use MFD_Aux_module
  use MFD_module

  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"     

  PetscInt :: local_id, ghosted_cell_id

  PetscReal, pointer :: porosity_loc_p(:), &
                        volume_p(:), &
                        accum_p(:), xx_faces_loc_p(:)

  type(realization_type) :: realization


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal :: mass_conserv, res(1)
  PetscErrorCode :: ierr

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  PetscScalar :: sq_faces(6), den(6), dden_dp(6), max_violation
  PetscInt :: jface, j, numfaces, ghost_face_id, local_face_id, dir, bound_id
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  field => realization%field

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_loc_p, ierr)

     max_violation = 0.

     do local_id = 1, grid%nlmax 
!     do local_id = 1, 1
        if ((local_id.ne.1).and.(local_id.ne.150)) cycle


        mass_conserv = 0.

        aux_var => grid%MFD%aux_vars(local_id)
        ghosted_cell_id = grid%nL2G(local_id)
        do j = 1, aux_var%numfaces
           ghost_face_id = aux_var%face_id_gh(j)
           local_face_id = grid%fG2L(ghost_face_id)
           conn => grid%faces(ghost_face_id)%conn_set_ptr
           jface = grid%faces(ghost_face_id)%id
           sq_faces(j) = conn%area(jface)
           call MFDComputeDensity(global_aux_vars(ghosted_cell_id), xx_faces_loc_p(ghost_face_id), &
                                                      den(j), dden_dp(j), option) 

#ifdef DASVYAT_DEBUG
           write(*,*)  local_face_id, jface, conn%itype, & 
                         conn%cntr(1,jface), conn%cntr(2,jface), conn%cntr(3,jface), xx_faces_loc_p(ghost_face_id)
#endif

           dir = 1

           if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(jface)==ghosted_cell_id) dir = -1 


           if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
               if (local_face_id > 0) then
                   bound_id = grid%fL2B(local_face_id)
                   if (bound_id>0) then
                       mass_conserv = mass_conserv - den(j)*patch%boundary_velocities(option%nphase, bound_id)*sq_faces(j)
#ifdef DASVYAT_DEBUG
                       write(*,*) "flux bnd", -patch%boundary_velocities(option%nphase, bound_id),&
                                      - den(j)*patch%boundary_velocities(option%nphase, bound_id)*sq_faces(j)
#endif
                   end if
               end if
           else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
               mass_conserv = mass_conserv + den(j)*patch%internal_velocities(option%nphase, jface)*sq_faces(j)* dir
#ifdef DASVYAT_DEBUG
                       write(*,*) "flux int",  patch%internal_velocities(option%nphase, jface)* dir, &
                                   den(j)*patch%internal_velocities(option%nphase, jface)*sq_faces(j)* dir
#endif
           end if
         end do
        
         call RichardsAccumulation(rich_aux_vars(ghosted_cell_id), &
                                global_aux_vars(ghosted_cell_id), &
                                porosity_loc_p(ghosted_cell_id), &
                                volume_p(local_id), &
                                option, res)

         mass_conserv = mass_conserv + res(1) - accum_p(local_id)
   
         max_violation = max(max_violation, abs(mass_conserv))

#ifdef DASVYAT_DEBUG
         write(*,*) "accumulation ", res(1) - accum_p(local_id)

         write(*,*) "Mass conservation ", local_id, mass_conserv
#endif
         
     end do

#ifdef DASVYAT_DEBUG
     write(*,*) "Mass conservation violation", max_violation
#endif
 
     call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
     call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
     call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
     call VecRestoreArrayF90(field%flow_accum, xx_faces_loc_p, ierr)

!     stop
#ifdef DASVYAT_DEBUG
     read(*,*)
#endif

end subroutine RichardsCheckMassBalancePatch

! ************************************************************************** !
!
! RichardsZeroMassBalDeltaPatch: Zeros mass balance delta array
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RichardsZeroMassBalDeltaPatch(realization)
 
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
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_aux_vars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine RichardsZeroMassBalDeltaPatch

! ************************************************************************** !
!
! RichardsUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RichardsUpdateMassBalancePatch(realization)
 
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
  do iconn = 1, patch%aux%Richards%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance = &
      patch%aux%Global%aux_vars(iconn)%mass_balance + &
      patch%aux%Global%aux_vars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Richards%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance = &
        global_aux_vars_bc(iconn)%mass_balance + &
        global_aux_vars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

  if (patch%aux%Richards%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Richards%num_aux_ss
      global_aux_vars_ss(iconn)%mass_balance = &
        global_aux_vars_ss(iconn)%mass_balance + &
        global_aux_vars_ss(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif

end subroutine RichardsUpdateMassBalancePatch

! ************************************************************************** !
!
! RichardsUpdateAuxVars: Updates the auxiliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateAuxVars(realization)

  use Realization_class
  use Grid_module, only : STRUCTURED_GRID_MIMETIC

  type(realization_type) :: realization
  
  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
!          if (.not.cur_patch%aux%Richards%aux_vars_cell_pressures_up_to_date) then
!             call RichardsUpdateCellPressurePatch(realization)
!          end if
    call RichardsUpdateAuxVarsPatchMFDLP(realization)
  else
    call RichardsUpdateAuxVarsPatch(realization)
  end if  

end subroutine RichardsUpdateAuxVars


! ************************************************************************** !
!
! RichardsInitialPressureReconstruction: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



! subroutine RichardsInitialPressureReconstruction(realization)
! 
!   use Realization_class
!   use Patch_module
! 
!   type(realization_type) :: realization
!   
!   type(level_type), pointer :: cur_level
!   type(patch_type), pointer :: cur_patch
! 
!   
!   cur_level => realization%level_list%first
!   do
!     if (.not.associated(cur_level)) exit
!     cur_patch => cur_level%patch_list%first
!     do
!       if (.not.associated(cur_patch)) exit
!       realization%patch => cur_patch
!       call RichardsInitialPressureReconstructionPatch(realization)
!       cur_patch => cur_patch%next
!     enddo
!     cur_level => cur_level%next
!   enddo
! 
! end subroutine RichardsInitialPressureReconstruction



! ************************************************************************** !
!
! RichardsUpdateCellPressure: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



subroutine RichardsUpdateCellPressure(realization)

  use Realization_class

  type(realization_type) :: realization

  if (.not.realization%patch%aux%Richards% &
      aux_vars_cell_pressures_up_to_date) then
    call RichardsUpdateCellPressurePatch(realization)
  end if

end subroutine RichardsUpdateCellPressure 

! ************************************************************************** !
!
! RichardsUpdateAuxVarsPatch: Updates the auxiliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateAuxVarsPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:) 
  type(richards_auxvar_type), pointer :: rich_aux_vars_bc(:)
  type(richards_auxvar_type), pointer :: rich_aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)  
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)  
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase, i
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)  
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  Vec :: phi
  
  call PetscLogEventBegin(logging%event_r_auxvars,ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  rich_aux_vars_ss => patch%aux%Richards%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
    
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    call RichardsAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id),rich_aux_vars(ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
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
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC)
          xxbc(1) = xx_loc_p(ghosted_id)
      end select
     
 
      call RichardsAuxVarCompute(xxbc(1),rich_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
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

      call RichardsAuxVarCopy(rich_aux_vars(ghosted_id), &
                              rich_aux_vars_ss(sum_connection),option)
      call GlobalAuxVarCopy(global_aux_vars(ghosted_id), &
                            global_aux_vars_ss(sum_connection),option)

    enddo
    source_sink => source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  

  ! Compute gradient using a least squares approach at each control volume
  if(realization%discretization%lsm_flux_method) then
    call RichardsUpdateLSMAuxVarsPatch(realization)
  endif

  patch%aux%Richards%aux_vars_up_to_date = PETSC_TRUE

  call PetscLogEventEnd(logging%event_r_auxvars_bc,ierr)

end subroutine RichardsUpdateAuxVarsPatch

! ************************************************************************** !
!> This routine computes the gradient of pressure using a least-square-method
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 11/24/12
! ************************************************************************** !
subroutine RichardsUpdateLSMAuxVarsPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:) 
  type(richards_auxvar_type), pointer :: rich_aux_vars_bc(:)
  type(richards_auxvar_type), pointer :: rich_aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)  
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)  
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase, i
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)  
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  PetscInt :: max_stencil_width
  PetscInt :: nid
  Vec :: phi,A,B
  PetscReal, pointer :: phi_p(:), b_p(:)
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal :: distance_gravity,distance_gravity0, den_avg
  PetscReal :: dx,dy,dz
  PetscReal :: P,Pn

  max_stencil_width = 2

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  rich_aux_vars_ss => patch%aux%Richards%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
    
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  

  select case(grid%itype)
    case(STRUCTURED_GRID)
      cell_neighbors => grid%structured_grid%cell_neighbors
    case(UNSTRUCTURED_GRID)
      option%io_buffer='GridComputeMinv() not implemented for unstructured grid.'
      call printErrMsg(option)
  end select

  do ghosted_id=1,grid%ngmax
    if(grid%ghosted_level(ghosted_id)<max_stencil_width) then
      call VecCreateSeq(PETSC_COMM_SELF,cell_neighbors(0,ghosted_id),phi,ierr)
      call VecGetArrayF90(phi,phi_p,ierr)
      do nid=1,cell_neighbors(0,ghosted_id)
        distance_gravity = option%gravity(1)*grid%x(cell_neighbors(nid,ghosted_id)) + &
                           option%gravity(2)*grid%y(cell_neighbors(nid,ghosted_id)) + &
                           option%gravity(3)*grid%z(cell_neighbors(nid,ghosted_id))

        Pn = xx_loc_p(cell_neighbors(nid,ghosted_id)) - &
             global_aux_vars(cell_neighbors(nid,ghosted_id))%den(1)*FMWH2O*distance_gravity

        distance_gravity0 = option%gravity(1)*grid%x(ghosted_id) + &
                            option%gravity(2)*grid%y(ghosted_id) + &
                            option%gravity(3)*grid%z(ghosted_id)
        P = xx_loc_p(ghosted_id) - &
             global_aux_vars(ghosted_id)%den(1)*FMWH2O*distance_gravity0

        phi_p(nid) = Pn-P
      enddo

      do nid=1,cell_neighbors(0,ghosted_id)
        if(abs(phi_p(nid))<1.D-20) phi_p(nid) = 0.d0
      enddo

      call VecRestoreArrayF90(phi,phi_p,ierr)

      call VecCreateSeq(PETSC_COMM_SELF,3,A,ierr)
      call MatMult(grid%dispT(ghosted_id),phi,A,ierr)

      call VecCreateSeq(PETSC_COMM_SELF,3,B,ierr)
      call MatMult(grid%Minv(ghosted_id),A,B,ierr)
      call VecGetArrayF90(B,b_p,ierr)

      ! If the gradient is too small, make it zero
      if(abs(b_p(1))<1.D-20) b_p(1) = 0.d0
      if(abs(b_p(2))<1.D-20) b_p(2) = 0.d0
      if(abs(b_p(3))<1.D-20) b_p(3) = 0.d0

      global_aux_vars(ghosted_id)%dphi(1,:) = b_p(:)

      call VecRestoreArrayF90(B,b_p,ierr)

      call VecDestroy(phi,ierr)
      call VecDestroy(A,ierr)
      call VecDestroy(B,ierr)
    endif
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

end subroutine RichardsUpdateLSMAuxVarsPatch

! ! ************************************************************************** !
! !
! ! RichardsInitialPressureReconstructionPatch: Computes cell-centered pressures and 
! ! author: Daniil Svyatskiy
! ! date: 12/01/10
! !
! ! ************************************************************************** !
! 
! 
! subroutine RichardsInitialPressureReconstructionPatch(realization)
! 
! 
!   use Realization_class
!   use Discretization_module
!   use Patch_module
!   use Option_module
!   use Field_module
!   use Grid_module
!   use Connection_module
!   use MFD_module
!   use MFD_aux_module
!   
!   implicit none
! 
!   type(realization_type) :: realization
!   type(discretization_type), pointer :: discretization
! 
!   
!   type(option_type), pointer :: option
!   type(patch_type), pointer :: patch
!   type(grid_type), pointer :: grid
!   type(field_type), pointer :: field
!   type(connection_set_type), pointer :: cur_connection_set
!   type(richards_auxvar_type), pointer :: rich_aux_vars(:)
!   type(global_auxvar_type), pointer :: global_aux_vars(:)
!   type(mfd_auxvar_type), pointer :: mfd_aux_var
! 
!   PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
!   PetscInt :: iphasebc, iphase
!   PetscInt :: ghost_face_id, iface, jface, numfaces
!   PetscReal, pointer :: xx_p(:), xx_loc_faces_p(:)
!   PetscReal, pointer :: sq_faces(:), faces_pr(:)
!   PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
!   PetscErrorCode :: ierr
! 
! 
! 
! #ifdef DASVYAT
!   write(*,*) "ENTER RichardsInitialPressureReconstructionPatch"
! #endif
! 
! 
!   discretization => realization%discretization
!   option => realization%option
!   patch => realization%patch
!   grid => patch%grid
!   field => realization%field
! 
! 
!   rich_aux_vars => patch%aux%Richards%aux_vars
!   global_aux_vars => patch%aux%Global%aux_vars
! 
!   call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces, ierr)
! 
!     
!   call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
!   call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
! 
! 
!   numfaces = 6     ! hex only
!   allocate(sq_faces(numfaces))
!   allocate(faces_pr(numfaces))
! 
! 
!   do local_id = 1, grid%nlmax
! 
!     ghosted_id = grid%nL2G(local_id)   
!  
!     !geh - Ignore inactive cells with inactive materials
!     if (patch%imat(ghosted_id) <= 0) cycle
!     mfd_aux_var => grid%MFD%aux_vars(local_id)
!     do j = 1, mfd_aux_var%numfaces
!        ghost_face_id = mfd_aux_var%face_id_gh(j)
!        cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
!        jface = grid%faces(ghost_face_id)%id
!        sq_faces(j) = cur_connection_set%area(jface)
! 
! !       if (cur_connection_set%itype == INTERNAL_CONNECTION_TYPE) then
! !         xx_loc_faces_p(ghost_face_id) = 1. !  DASVYAT test
! !       end if 
! 
!        faces_pr(j) = xx_loc_faces_p(ghost_face_id)
!     end do
!  
!         
!         !geh - Ignore inactive cells with inactive materials
!         if (patch%imat(ghosted_id) <= 0) cycle
!  
!         Res(1) = 0
! 
!         source_f = 0.
! !        Res(1) = 0.
! 
!         call MFDAuxReconstruct(faces_pr, source_f, mfd_aux_var,&
!                                rich_aux_vars(ghosted_id),global_aux_vars(ghosted_id), Res, &
!                                sq_faces, option, xx_p(local_id:local_id))
! 
!    
! !       write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)
! 
!   enddo
! 
! 
!  patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_TRUE 
! 
! 
!   deallocate(sq_faces)
!   deallocate(faces_pr)
! 
! 
!   call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
!   call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
! 
!   call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF) 
! 
! 
! !  call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_xx_loc_faces, field%flow_xx_faces, &   ! DASVYAT test
! !                                INSERT_VALUES,SCATTER_REVERSE, ierr)
! !   call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_xx_loc_faces, field%flow_xx_faces,&    ! DASVYAT test 
! !                                INSERT_VALUES,SCATTER_REVERSE, ierr)
! 
! 
!   
!   write(*,*) "DiscretizationGlobalToLocal"
! !   read(*,*)
!     stop
! 
! end subroutine RichardsInitialPressureReconstructionPatch

! ************************************************************************** !
!
! RichardsUpdateCellPressurePatch: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



subroutine RichardsUpdateCellPressurePatch(realization)


  use Realization_class
  use Discretization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use MFD_module
  use MFD_aux_module
  
  implicit none

  type(realization_type) :: realization

  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(mfd_auxvar_type), pointer :: mfd_aux_var

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
  PetscInt :: iphasebc, iphase
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_p(:), xx_loc_faces_p(:), work_loc_faces_p(:)
  PetscReal, pointer ::  volume_p(:), porosity_p(:), accum_p(:)
  PetscReal, pointer :: sq_faces(:), faces_DELTA_pr(:), faces_pr(:) 
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr



#ifdef DASVYAT_DEBUG
  write(*,*) "ENTER RichardsUpdateCellPressurePatch"
#endif


  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  discretization => realization%discretization


  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

    
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecGetArrayF90(field%work_loc_faces, work_loc_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity0, porosity_p,ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p,ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p,ierr)


  numfaces = 6     ! hex only
  allocate(sq_faces(numfaces))
  allocate(faces_DELTA_pr(numfaces))
  allocate(faces_pr(numfaces))


  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)   
 
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    mfd_aux_var => grid%MFD%aux_vars(local_id)
    do j = 1, mfd_aux_var%numfaces
       ghost_face_id = mfd_aux_var%face_id_gh(j)
       cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
       jface = grid%faces(ghost_face_id)%id
       sq_faces(j) = cur_connection_set%area(jface)
       faces_pr(j) = work_loc_faces_p(ghost_face_id)
       
       faces_pr(j) = xx_loc_faces_p(ghost_face_id) 
       
       faces_DELTA_pr(j) = xx_loc_faces_p(ghost_face_id) - work_loc_faces_p(ghost_face_id)
!       faces_DELTA_pr(j) = 0
    end do
 
        
        !geh - Ignore inactive cells with inactive materials
        if (patch%imat(ghosted_id) <= 0) cycle
        

       call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_p(local_id), &
                                volume_p(local_id), &
                                option,Res)


        Res(1) = Res(1) - accum_p(local_id)

        source_f = 0.
!        Res(1) = 0.
!       write(*,*) "Before pres", xx_p(local_id)

        call MFDAuxUpdateCellPressure(faces_pr, faces_DELTA_pr, mfd_aux_var,&
                               option, xx_p(local_id:local_id))

   
!       write(*,*) "After pres", xx_p(local_id)

  enddo

 !     write(*,*) "RichardsUpdateCellPressurePatch"
 !     read(*,*)

 patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_TRUE 


  deallocate(sq_faces)
  deallocate(faces_DELTA_pr)
  deallocate(faces_pr)

  call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF)


  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecRestoreArrayF90(field%work_loc_faces, work_loc_faces_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity0, porosity_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p,ierr)


end subroutine RichardsUpdateCellPressurePatch




! ************************************************************************** !
!
! RichardsUpdateAuxVarsPatchMFDLP: Computes  updates the auxiliary variables associated with 
!                        the Richards problem for LP formulation
! author: Daniil Svyatskiy
! date: 07/29/10
!
! ************************************************************************** !

subroutine RichardsUpdateAuxVarsPatchMFDLP(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  use MFD_module
  use MFD_aux_module
  
  implicit none

  type(realization_type) :: realization

#ifdef DASVYAT
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(mfd_auxvar_type), pointer :: mfd_aux_var

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
  PetscInt :: iphasebc, iphase, LP_cell_id, LP_face_id
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_p(:), xx_LP_loc_p(:), bc_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)  
  PetscReal, pointer :: sq_faces(:), darcy_v(:), faces_pr(:)
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr


! For Boundary Conditions
  PetscReal :: xxbc(realization%option%nflowdof)
  type(richards_auxvar_type), pointer :: rich_aux_vars_bc(:) 
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)


#ifdef DASVYAT
!   write(*,*) "ENTER RichardsUpdateAuxVarsPatchMFDLP"
#endif

  call PetscLogEventBegin(logging%event_r_auxvars,ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field


  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

    
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_LP_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  




  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

     !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

 
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    LP_cell_id = grid%ngmax_faces + ghosted_id

!     write(*,*) "press", "xx", xx_loc_p(local_id), "xx_LP", xx_LP_p(LP_cell_id)
   
       call RichardsAuxVarCompute(xx_LP_loc_p(LP_cell_id:LP_cell_id),rich_aux_vars(ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)

       local_id = grid%nG2L(ghosted_id)
       if (local_id > 0) xx_p(local_id) = xx_LP_loc_p(LP_cell_id)

!        write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)

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
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          xxbc(1) = xx_LP_loc_p(grid%ngmax_faces + ghosted_id)
        case(UNIT_GRADIENT_BC)
          ! the auxvar is not needed for unit gradient
          cycle
      end select
      
      call RichardsAuxVarCompute(xxbc(1),rich_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call PetscLogEventEnd(logging%event_r_auxvars,ierr)

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_LP_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p,ierr)  

#if 0
  write(*,*) "EXIT RichardsUpdateAuxVarsPatchMFDLP"
  read(*,*)
#endif

#endif

end subroutine RichardsUpdateAuxVarsPatchMFDLP

! ************************************************************************** !
!
! RichardsInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RichardsInitializeTimestep(realization)

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
!
! RichardsUpdateSolution: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine RichardsUpdateSolution(realization)

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
!
! RichardsUpdateSolutionPatch: Updates data in module after a successful time 
!                             step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine RichardsUpdateSolutionPatch(realization)

  use Realization_class
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call RichardsUpdateMassBalancePatch(realization)
  endif

end subroutine RichardsUpdateSolutionPatch

! ************************************************************************** !
!
! RichardsUpdateFixedAccum: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateFixedAccum(realization)

  use Realization_class

  type(realization_type) :: realization
  
  call RichardsUpdateFixedAccumPatch(realization)

end subroutine RichardsUpdateFixedAccum

! ************************************************************************** !
!
! RichardsUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                accumulation term
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateFixedAccumPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use MFD_aux_module
  use MFD_module
  use Connection_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)

  PetscInt :: ghosted_id, local_id, numfaces, jface, ghost_face_id, j
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          accum_p(:), perm_xx_loc_p(:)
!  PetscReal, pointer :: xx_faces_p(:)

!  PetscReal, pointer :: faces_pr(:), sq_faces(:)

!  type(mfd_auxvar_type), pointer :: mfd_aux_var
!  type(connection_set_type), pointer :: cur_connection_set
!  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)                        
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
    
  call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)  

!  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)

  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)

!  numfaces = 6     ! hex only
!  allocate(sq_faces(numfaces))
!  allocate(faces_pr(numfaces))

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsAuxVarCompute(xx_p(local_id:local_id), &
                   rich_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                   patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                   porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                        
                   option)
    call RichardsAccumulation(rich_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option,accum_p(local_id:local_id))
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tortuosity_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)  


!  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)

  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)

#if 0
!  call RichardsNumericalJacTest(field%flow_xx,realization)
#endif

end subroutine RichardsUpdateFixedAccumPatch

! ************************************************************************** !
!
! RichardsNumericalJacTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsNumericalJacTest(xx,realization)

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
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
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
      call GridVecRestoreArrayF90(grid,res_pert,vec_p,ierr)
!    enddo
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

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
!
! RichardsAccumDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsAccumDerivative(rich_aux_var,global_aux_var,por,vol, &
                                   option,sat_func,J)

  use Option_module
  use Saturation_Function_module
  
  implicit none

  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: vol, por
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: vol_over_dt
  PetscReal :: tempreal 

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert
  type(global_auxvar_type) :: global_aux_var_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  vol_over_dt = vol/option%flow_dt
      
!#define USE_COMPRESSIBLITY
#ifndef USE_COMPRESSIBLITY  
  ! accumulation term units = dkmol/dp
  J(1,1) = (global_aux_var%sat(1)*rich_aux_var%dden_dp+ &
            rich_aux_var%dsat_dp*global_aux_var%den(1))* &
           por*vol_over_dt

#else
  tempreal = exp(-1.d-7*(global_aux_var%pres(1)-option%reference_pressure))
  J(1,1) = ((global_aux_var%sat(1)*rich_aux_var%dden_dp+ &
             rich_aux_var%dsat_dp*global_aux_var%den(1))* &
            1.d0-(1.d0-por)*tempreal + &
            global_aux_var%sat(1)*global_aux_var%den(1)* &
            (por-1.d0)*-1.d-7*tempreal)* &
           vol_over_dt
#endif  
  
  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert,option)  
    call RichardsAuxVarCopy(rich_aux_var,rich_aux_var_pert,option)
    call GlobalAuxVarCopy(global_aux_var,global_aux_var_pert,option)
    x(1) = global_aux_var%pres(1)
    call RichardsAccumulation(rich_aux_var,global_aux_var,por,vol,option,res)
    ideriv = 1
    pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
    x_pert = x
    if (x_pert(ideriv) < option%reference_pressure) pert = -1.d0*pert
    x_pert(ideriv) = x_pert(ideriv) + pert
    
    call RichardsAuxVarCompute(x_pert(1),rich_aux_var_pert,global_aux_var_pert, &
                               sat_func,0.d0,0.d0,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', aux_var%dvis_dp, (aux_var_pert%vis-aux_var%vis)/pert(ideriv)
!         print *, 'dkr_dp:', aux_var%dkr_dp, (aux_var_pert%kr-aux_var%kr)/pert(ideriv)
          print *, 'dsat_dp:', aux_var%dsat_dp, (global_aux_var_pert%sat-global_aux_var%sat)/pert
          print *, 'dden_dp:', aux_var%dden_dp, (global_aux_var_pert%den-global_aux_var%den)/pert
!          print *, 'dkvr_dp:', aux_var%dkvr_dp, (rich_aux_var_pert%kvr-rich_aux_var%kvr)/pert
          print *, 'dkvr_x_dp:', aux_var%dkvr_x_dp, (rich_aux_var_pert%kvr_x-rich_aux_var%kvr_x)/pert
          print *, 'dkvr_y_dp:', aux_var%dkvr_y_dp, (rich_aux_var_pert%kvr_y-rich_aux_var%kvr_y)/pert
          print *, 'dkvr_z_dp:', aux_var%dkvr_z_dp, (rich_aux_var_pert%kvr_z-rich_aux_var%kvr_z)/pert
      end select     
#endif     
    call RichardsAccumulation(rich_aux_var_pert,global_aux_var_pert,por,vol, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarStrip(global_aux_var_pert)  
  endif
   
end subroutine RichardsAccumDerivative


! ************************************************************************** !
!
! RichardsAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsAccumulation(rich_aux_var,global_aux_var,por,vol, &
                                option,Res)

  use Option_module
  
  implicit none

  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: vol, por, por1
       
  ! accumulation term units = kmol/s

#ifdef DASVYAT
!  write(*,*) "sat ",global_aux_var%sat(1), " den ", global_aux_var%den(1)
!  write(*,*) "por ", por, " volume ", vol, " dt ", option%flow_dt
#endif

#ifndef USE_COMPRESSIBILITY
    Res(1) = global_aux_var%sat(1) * global_aux_var%den(1) * por * vol / &
           option%flow_dt
#else
    por1 = 1.d0-(1.d0-por)*exp(-1.d-7*(global_aux_var%pres(1)- &
                                       option%reference_pressure))
    Res(1) = global_aux_var%sat(1) * global_aux_var%den(1) * por1 * vol / &
           option%flow_dt
#endif
    
end subroutine RichardsAccumulation

! ************************************************************************** !
!
! RichardsFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsFluxDerivative(rich_aux_var_up,global_aux_var_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  option,sat_func_up,sat_func_dn,Jup,Jdn)
  use Option_module 
  use Saturation_Function_module                        
  
  implicit none
  
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy, area, dist(3)
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
  PetscReal :: q
  PetscReal :: ukvr,Dq
!  PetscReal :: ukvr_x, ukvr_y, ukvr_z, Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
!  PetscReal :: dukvr_x_dp_up, dukvr_x_dp_dn
!  PetscReal :: dukvr_y_dp_up, dukvr_y_dp_dn
!  PetscReal :: dukvr_z_dp_up, dukvr_z_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert_up, rich_aux_var_pert_dn
  type(global_auxvar_type) :: global_aux_var_pert_up, global_aux_var_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
  v_darcy = 0.D0 
  ukvr = 0.d0
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)
    dden_ave_dp_up = upweight*rich_aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_aux_var_dn%dden_dp

    gravity = (upweight*global_aux_var_up%den(1) + &
              (1.D0-upweight)*global_aux_var_dn%den(1)) &
              * FMWH2O * dist_gravity
    dgravity_dden_up = upweight*FMWH2O*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*rich_aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_aux_var_dn%dden_dp

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dist(1))==1) then
        ukvr = rich_aux_var_up%kvr_x
        dukvr_dp_up = rich_aux_var_up%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_aux_var_up%kvr_y
        dukvr_dp_up = rich_aux_var_up%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_aux_var_up%kvr_z
        dukvr_dp_up = rich_aux_var_up%dkvr_z_dp
      end if
#else
      ukvr = rich_aux_var_up%kvr
      dukvr_dp_up = rich_aux_var_up%dkvr_dp
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY    
      if (dabs(dist(1))==1) then
        ukvr = rich_aux_var_dn%kvr_x
        dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_aux_var_dn%kvr_y
        dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_aux_var_dn%kvr_z
        dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
      end if
#else
      ukvr = rich_aux_var_dn%kvr
      dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
    endif
  endif

 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call RichardsAuxVarCopy(rich_aux_var_up,rich_aux_var_pert_up,option)
    call RichardsAuxVarCopy(rich_aux_var_dn,rich_aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    x_up(1) = global_aux_var_up%pres(1)
    x_dn(1) = global_aux_var_dn%pres(1)
    call RichardsFlux(rich_aux_var_up,global_aux_var_up,por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_dn,global_aux_var_dn,por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res)
    ideriv = 1
!    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_up = max(dabs(x_up(ideriv)*perturbation_tolerance),0.1d0)
    if (x_up(ideriv) < option%reference_pressure) pert_up = -1.d0*pert_up
!    pert_dn = x_dn(ideriv)*perturbation_tolerance
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_up = x_up
    x_pert_dn = x_dn
    x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    call RichardsAuxVarCompute(x_pert_up(1),rich_aux_var_pert_up, &
                               global_aux_var_pert_up,sat_func_up, &
                               0.d0,0.d0,option)
    call RichardsAuxVarCompute(x_pert_dn(1),rich_aux_var_pert_dn, &
                               global_aux_var_pert_dn,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsFlux(rich_aux_var_pert_up,global_aux_var_pert_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_dn,global_aux_var_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res_pert_up)
    call RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_pert_dn,global_aux_var_pert_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)    
  endif

end subroutine RichardsFluxDerivative

! ************************************************************************** !
!> This routine computes the derivatives of the internal flux terms for the
!! Jacobian, when least-squares-method is used to compute the gradients.
!!
!! NOTE: Each internal flux term, depends on states of 'up', 'down', and
!!       neighbors of 'up' & 'down' cells. Thus, jacobian matrix should have
!!       entries in multiple columns of the Jacobian matrix that corresponds
!!       to neigboring cells for a given internal flux term. But, the entries
!!       corresponding to neigboring cells are neglected presently.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/17/12
! ************************************************************************** !
subroutine RichardsLSMFluxDerivative(rich_aux_var_up,global_aux_var_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  distance, &
                                  option,sat_func_up,sat_func_dn, &
                                  jacfac,ghosted_id_up,ghosted_id_dn, &
                                  cell_neighbors, &
                                  x,y,z, &
                                  bnd_cell, &
                                  Jup,Jdn)
  use Option_module
  use Saturation_Function_module

  implicit none

  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy, area, dist(3)
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: distance
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscReal,pointer :: jacfac(:,:,:)
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal, pointer :: x(:),y(:),z(:)
  PetscBool, pointer :: bnd_cell(:)
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: dphi_x_dp_up,dphi_y_dp_up,dphi_z_dp_up
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn

  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn
  PetscReal :: dHup_x_dp_up, dHup_y_dp_up, dHup_z_dp_up
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  PetscInt :: nid_up, nid_dn,ii
  type(richards_auxvar_type) :: rich_aux_var_pert_up, rich_aux_var_pert_dn
  type(global_auxvar_type) :: global_aux_var_pert_up, global_aux_var_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  v_darcy = 0.D0
  ukvr = 0.d0

  Jup = 0.d0
  Jdn = 0.d0

  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0

  nid_up = -1
  nid_dn = -1

  ! One of the cells is a boundary cell, so use 2-point flux
  if(bnd_cell(ghosted_id_up).or.bnd_cell(ghosted_id_dn)) then
    call RichardsFluxDerivative(rich_aux_var_up,global_aux_var_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  option,sat_func_up,sat_func_dn,Jup,Jdn)
    return
  endif

  do ii = 1,cell_neighbors(0,ghosted_id_dn)
    if (cell_neighbors(ii,ghosted_id_dn) == ghosted_id_up) nid_up = ii
  enddo
  do ii = 1,cell_neighbors(0,ghosted_id_up)
    if (cell_neighbors(ii,ghosted_id_up) == ghosted_id_dn) nid_dn = ii
  enddo
  if(nid_up<0 .or. nid_dn <0 ) then
    option%io_buffer = 'Neighbors not found'
    call printErrMsg(option)
  endif

  ! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then
      upweight=1.d0
    endif

    density_ave =        upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)
    dden_ave_dp_up = upweight*rich_aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_aux_var_dn%dden_dp

    dphi_x =       upweight *global_aux_var_up%dphi(1,1) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,1)
    dphi_y =       upweight *global_aux_var_up%dphi(1,2) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,2)
    dphi_z =       upweight *global_aux_var_up%dphi(1,3) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,3)

    dphi = -(dist(1)*dphi_x + dist(2)*dphi_y + dist(3)*dphi_z)*distance

    ! H_x = P - W * rho * g *z
    dHup_x_dp_up = 1.d0 - FMWH2O*rich_aux_var_up%dden_dp*option%gravity(1)*x(ghosted_id_up)
    dHup_y_dp_up = 1.d0 - FMWH2O*rich_aux_var_up%dden_dp*option%gravity(2)*y(ghosted_id_up)
    dHup_z_dp_up = 1.d0 - FMWH2O*rich_aux_var_up%dden_dp*option%gravity(3)*z(ghosted_id_up)

    dHdn_x_dp_dn = 1.d0 - FMWH2O*rich_aux_var_dn%dden_dp*option%gravity(1)*x(ghosted_id_dn)
    dHdn_y_dp_dn = 1.d0 - FMWH2O*rich_aux_var_dn%dden_dp*option%gravity(2)*y(ghosted_id_dn)
    dHdn_z_dp_dn = 1.d0 - FMWH2O*rich_aux_var_dn%dden_dp*option%gravity(3)*z(ghosted_id_dn)

    dphi_x_dp_up =       upweight *jacfac(ghosted_id_up,0     ,1)*dHup_x_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,1)*dHup_x_dp_up
    dphi_y_dp_up =       upweight *jacfac(ghosted_id_up,0     ,2)*dHup_y_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,2)*dHup_y_dp_up
    dphi_z_dp_up =       upweight *jacfac(ghosted_id_up,0     ,3)*dHup_z_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,3)*dHup_z_dp_up

    dphi_x_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,1)*dHdn_x_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,1)*dHdn_x_dp_dn
    dphi_y_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,2)*dHdn_y_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,2)*dHdn_y_dp_dn
    dphi_z_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,3)*dHdn_z_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,3)*dHdn_z_dp_dn

    dphi_dp_up = -(dist(1)*dphi_x_dp_up + dist(2)*dphi_y_dp_up + dist(3)*dphi_z_dp_up)
    dphi_dp_dn = -(dist(1)*dphi_x_dp_dn + dist(2)*dphi_y_dp_dn + dist(3)*dphi_z_dp_dn)

    if (dphi>=0.D0) then
      ukvr = rich_aux_var_up%kvr
      dukvr_dp_up = rich_aux_var_up%dkvr_dp
    else
      ukvr = rich_aux_var_dn%kvr
      dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
    endif

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi

      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area

      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

    endif

  endif

end subroutine RichardsLSMFluxDerivative

! ************************************************************************** !
!
! RichardsFlux: Computes the internal flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                        por_up,sir_up,dd_up,perm_up, &
                        rich_aux_var_dn,global_aux_var_dn, &
                        por_dn,sir_dn,dd_dn,perm_dn, &
                        area, dist, dist_gravity,upweight, &
                        option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy,area, dist(3)
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
  fluxm = 0.d0
  v_darcy = 0.D0  
  ukvr = 0.d0
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    

    density_ave = upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)

    gravity = (upweight*global_aux_var_up%den(1) + &
              (1.D0-upweight)*global_aux_var_dn%den(1)) &
              * FMWH2O * dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)  + gravity

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_z
      end if
#else
      ukvr = rich_aux_var_up%kvr
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(norm(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
#else
      ukvr = rich_aux_var_dn%kvr
#endif
    endif      

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area

      fluxm = q*density_ave       
    endif
  endif 

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsFlux

! ************************************************************************** !
!> This routine computes the internal flux terms for the residual term using
!! a least-square-method for gradient computation.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 12/10/12
! ************************************************************************** !
subroutine RichardsLSMFlux(rich_aux_var_up,global_aux_var_up, &
                           por_up,sir_up,dd_up,perm_up, &
                           rich_aux_var_dn,global_aux_var_dn, &
                           por_dn,sir_dn,dd_dn,perm_dn, &
                           area, dist, dist_gravity,upweight, distance, &
                           option, &
                           ghosted_id_up, ghosted_id_dn, cell_neighbors, &
                           bnd_cell, &
                           v_darcy,Res)
  use Option_module

  implicit none

  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy,area, dist(3)
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscInt  :: ghosted_id_up,ghosted_id_dn
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscBool, pointer :: bnd_cell(:)

  PetscInt :: ispec
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,dphi
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: fraction_upwind,distance

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.d0
  v_darcy = 0.D0
  ukvr = 0.d0

 ! One of the cells is a boundary cell, so use 2-point flux
  if(bnd_cell(ghosted_id_up).or.bnd_cell(ghosted_id_dn)) then
    call RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_dn,global_aux_var_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,Res)
    return
  endif

  ! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then
      upweight=1.d0
    endif

    density_ave = upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)

    dphi_x =       upweight *global_aux_var_up%dphi(1,1) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,1)
    dphi_y =       upweight *global_aux_var_up%dphi(1,2) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,2)
    dphi_z =       upweight *global_aux_var_up%dphi(1,3) + &
             (1.d0-upweight)*global_aux_var_dn%dphi(1,3)

    dphi = -(dist(1)*dphi_x + dist(2)*dphi_y + dist(3)*dphi_z)*distance

    if (dphi>=0.D0) then
      ukvr = rich_aux_var_up%kvr
    else
      ukvr = rich_aux_var_dn%kvr
    endif

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
      q = v_darcy * area
      fluxm = q*density_ave
    endif
  endif

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL

end subroutine RichardsLSMFlux

! ************************************************************************** !
!
! RichardsBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsBCFluxDerivative(ibndtype,aux_vars, &
                                    rich_aux_var_up,global_aux_var_up, &
                                    rich_aux_var_dn,global_aux_var_dn, &
                                    por_dn,sir_dn,perm_dn, &
                                    area, dist,option, &
                                    sat_func_dn,Jdn)
  use Option_module
  use Saturation_Function_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn
  PetscReal :: area
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(0:3)
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscReal :: v_darcy
  PetscReal :: q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: dden_ave_dp_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn
  PetscReal :: dukvr_dp_dn
  PetscReal :: dq_dp_dn
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert_dn, rich_aux_var_pert_up
  type(global_auxvar_type) :: global_aux_var_pert_dn, global_aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  
  v_darcy = 0.d0
  ukvr = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_dn = 0.d0
        
  ! Flow
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = aux_vars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_aux_var_up%den(1) + (1.D0-upweight)*global_aux_var_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*rich_aux_var_dn%dden_dp

        gravity = (upweight*global_aux_var_up%den(1) + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)) &
                  * FMWH2O * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_aux_var_dn%dden_dp

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY  
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_x
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_y
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_z
          end if
#else
          ukvr = rich_aux_var_up%kvr
#endif
        else
#ifdef USE_ANISOTROPIC_MOBILITY
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_x
            dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_y
            dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_z
            dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
          end if
#else
          ukvr = rich_aux_var_dn%kvr
          dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi + ukvr*dphi_dp_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
        else 
          density_ave = global_aux_var_dn%den(1)
          dden_ave_dp_dn = rich_aux_var_dn%dden_dp
        endif 
        q = v_darcy * area
      endif

    case(UNIT_GRADIENT_BC)

      if (global_aux_var_dn%sat(1) > sir_dn) then

        dphi = dot_product(option%gravity,dist(1:3))*global_aux_var_dn%den(1)*FMWH2O      
        density_ave = global_aux_var_dn%den(1)

        dphi_dp_dn = dot_product(option%gravity,dist(1:3))*rich_aux_var_dn%dden_dp*FMWH2O
        dden_ave_dp_dn = rich_aux_var_dn%dden_dp

        ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
        if (dabs(dabs(dist(1))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_x
          dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
        else if (dabs(dabs(dist(2))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_y
          dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
        else if (dabs(dabs(dist(3))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_z
          dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
        end if
#else
        ukvr = rich_aux_var_dn%kvr
        dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
     
        if (ukvr*perm_dn>floweps) then
          v_darcy = perm_dn * ukvr * dphi
          q = v_darcy*area
          dq_dp_dn = perm_dn*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif 
     
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call RichardsAuxVarCopy(rich_aux_var_up,rich_aux_var_pert_up,option)
    call RichardsAuxVarCopy(rich_aux_var_dn,rich_aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    x_up(1) = global_aux_var_up%pres(1)
    x_dn(1) = global_aux_var_dn%pres(1)
    ideriv = 1
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_up(ideriv) = x_dn(ideriv)
    endif
    call RichardsBCFlux(ibndtype,aux_vars, &
                        rich_aux_var_up,global_aux_var_up, &
                        rich_aux_var_dn,global_aux_var_dn, &
                        por_dn,sir_dn,perm_dn, &
                        area,dist,option,v_darcy,res)
    if (pressure_bc_type == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    ideriv = 1
!    pert_dn = x_dn(ideriv)*perturbation_tolerance    
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_dn = x_dn
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    x_pert_up = x_up
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_pert_up(ideriv) = x_pert_dn(ideriv)
    endif   
    call RichardsAuxVarCompute(x_pert_dn(1),rich_aux_var_pert_dn, &
                               global_aux_var_pert_dn,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsAuxVarCompute(x_pert_up(1),rich_aux_var_pert_up, &
                               global_aux_var_pert_up,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsBCFlux(ibndtype,aux_vars, &
                        rich_aux_var_pert_up,global_aux_var_pert_up, &
                        rich_aux_var_pert_dn,global_aux_var_pert_dn, &
                        por_dn,sir_dn,perm_dn, &
                        area,dist,option,v_darcy,res_pert_dn)
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)      
  endif

end subroutine RichardsBCFluxDerivative

! ************************************************************************** !
!
! RichardsBCFlux: Computes the  boundary flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsBCFlux(ibndtype,aux_vars, &
                          rich_aux_var_up, global_aux_var_up, &
                          rich_aux_var_dn, global_aux_var_dn, &
                          por_dn, sir_dn, perm_dn, &
                          area, dist, option,v_darcy,Res)
  use Option_module
#ifdef SURFACE_FLOW
  use water_eos_module
#endif
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn
  PetscReal :: v_darcy, area
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(0:3)
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm,q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x,dphi_y,dphi_z
#ifdef SURFACE_FLOW
  PetscReal :: rho, v_darcy_allowable
#endif
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  ukvr = 0.d0

  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = aux_vars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
   
        gravity = (upweight*global_aux_var_up%den(1) + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)) &
                  * FMWH2O * dist_gravity
       
        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
   
       if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_z
         end if
#else
         ukvr = rich_aux_var_up%kvr
#endif
       else
#ifdef USE_ANISOTROPIC_MOBILITY
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_z
         end if
#else
         ukvr = rich_aux_var_dn%kvr
#endif
       endif
        
       if (ukvr*Dq>floweps) then
        v_darcy = Dq * ukvr * dphi
#ifdef SURFACE_FLOW
        if(option%nsurfflowdof>0) then
          call density(option%reference_temperature,option%reference_pressure,rho)
          v_darcy_allowable = (global_aux_var_up%pres(1)-option%reference_pressure) &
                              /option%flow_dt/(-option%gravity(3))/rho
          v_darcy = min(v_darcy,v_darcy_allowable)
        endif
#endif
       endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)

        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
        else 
          density_ave = global_aux_var_dn%den(1)
        endif 
      endif

    case(UNIT_GRADIENT_BC)

      dphi = dot_product(option%gravity,dist(1:3))*global_aux_var_dn%den(1)*FMWH2O
      density_ave = global_aux_var_dn%den(1)

      ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(dist(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(dist(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
#else
      ukvr = rich_aux_var_dn%kvr
#endif
     
      if (ukvr*perm_dn>floweps) then
        v_darcy = perm_dn * ukvr * dphi
      endif 

  end select

  q = v_darcy * area

  fluxm = q*density_ave

  Res(1)=fluxm

end subroutine RichardsBCFlux

! ************************************************************************** !
!
! RichardsResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidual(snes,xx,r,realization,ierr)

  use Realization_class
  use Field_module
  use Discretization_module
  use Option_module
  use Logging_module

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
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  
  ! pass #1 for internal and boundary flux terms
  call RichardsResidualPatch1(snes,xx,r,realization,ierr)

  ! pass #2 for everything else
  call RichardsResidualPatch2(snes,xx,r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
!geh    call PetscViewerBinaryOpen(realization%option%mycomm,'Rresidual.bin',FILE_MODE_WRITE,&
!geh                              viewer,ierr)
!geh    call VecView(r,viewer,ierr)
!geh    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
!geh    call PetscViewerBinaryOpen(realization%option%mycomm,'Rxx.bin',FILE_MODE_WRITE,&
!geh                              viewer,ierr)
!geh    call VecView(xx,viewer,ierr)
!geh    call PetscViewerDestroy(viewer,ierr)

  endif

call PetscLogEventEnd(logging%event_r_residual,ierr)

#ifdef DASVYAT
!  write(*,*) "End of RichardsResidual"
!  read(*,*)
!  stop 
#endif
 
end subroutine RichardsResidual

! ************************************************************************** !
!
! RichardsResidualMFD: Computes the residual equation for MFD discretization
! author: Daniil Svyatskiy
! date: 05/26/10
!
! ************************************************************************** !
subroutine RichardsResidualMFD(snes,xx,r,realization,ierr)

  use Logging_module
  use Realization_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Connection_module

  implicit none

#include "finclude/petscsysdef.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
#ifdef DASVYAT

  PetscInt :: i, iface
  PetscReal, pointer :: flow_xx_loc_p(:), r_p(:)
  PetscReal :: rnorm
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: conn
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option


  call PetscLogEventBegin(logging%event_r_residual,ierr)

  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr)
  r_p = 0.
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr)

  call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r, &
                                INSERT_VALUES,SCATTER_REVERSE, ierr)
  call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r,&
                                INSERT_VALUES,SCATTER_REVERSE, ierr)


!   write(*,*) "Begin RichardsResidualMFD"
!   read(*,*)
 
  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
   call DiscretizationGlobalToLocalFaces(discretization, xx, field%flow_xx_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocalFaces(discretization, r, field%flow_r_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF)

!  write(*,*) "After DiscretizationGlobalToLocalFaces"
!  read(*,*)

  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)


  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xz_loc,field%perm_xz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xy_loc,field%perm_xy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yz_loc,field%perm_yz_loc,ONEDOF)


  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
         call RichardsUpdateCellPressure(realization)
  end if  

  
  ! pass #1 for flux terms and accumulation term
  call RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)

  ! pass #2 for boundary and source data
  call RichardsResidualPatchMFD2(snes,xx,r,realization,ierr)
!  call RichardsCheckMassBalancePatch(realization)

   call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces, ierr)

   call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r, &
                                ADD_VALUES,SCATTER_REVERSE, ierr)
   call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r,&
                                ADD_VALUES,SCATTER_REVERSE, ierr)

#if 0   
   call DiscretizationGlobalToLocalFaces(discretization, r, field%flow_r_loc_faces, NFLOWDOF)

  call GridVecGetArrayF90(realization%patch%grid, field%flow_r_loc_faces, r_p, ierr)

  call VecStrideNorm(r ,ZERO_INTEGER,NORM_INFINITY, rnorm ,ierr)


  call GridVecRestoreArrayF90(realization%patch%grid, field%flow_r_loc_faces, r_p, ierr)
#endif  

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RresidualMFD.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RxxMFD.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
 endif


  call PetscLogEventEnd(logging%event_r_residual,ierr)

#ifdef DASVYAT_DEBUG
   write(*,*) "End RichardsResidualMFD"
   read(*,*) 
#endif   
!   write(*,*) "End RichardsResidualMFD"
!   read(*,*) 
!   stop

#endif
  
end subroutine RichardsResidualMFD



subroutine RichardsResidualMFDLP(snes,xx,r,realization,ierr)

use Logging_module
  use Realization_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Connection_module

  implicit none

#include "finclude/petscsysdef.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr

#if DASVYAT

  PetscInt :: i, iface
  PetscReal, pointer :: flow_xx_loc_p(:), r_p(:)
  PetscReal :: rnorm
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: conn
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option


  call PetscLogEventBegin(logging%event_r_residual,ierr)

  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr)
  r_p = 0.
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr)

  call VecScatterBegin( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r, &
                                INSERT_VALUES,SCATTER_REVERSE, ierr)
  call VecScatterEnd ( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r,&
                                INSERT_VALUES,SCATTER_REVERSE, ierr)

!  if (option%test_res==3) then
!    if (realization%debug%vecview_residual) then
!      call PetscViewerASCIIOpen(realization%option%mycomm,'RresidualMFD_after.out', &
!                                viewer,ierr)
!      call VecView(r,viewer,ierr)
!      call PetscViewerDestroy(viewer,ierr)
!    endif
!    if (realization%debug%vecview_solution) then
!      call PetscViewerASCIIOpen(realization%option%mycomm,'RxxMFD_after.out', &
!                                viewer,ierr)
!      call VecView(xx,viewer,ierr)
!      call PetscViewerDestroy(viewer,ierr)
!   endif
!   write(*,*) "option%test_res==1"
!   do ierr = 1, 1000000
!   end do
 !  stop
!  end if

 
  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
   call DiscretizationGlobalToLocalLP(discretization, xx, field%flow_xx_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocalLP(discretization, r, field%flow_r_loc_faces, NFLOWDOF)
!    call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF)


!   call PetscViewerASCIIOpen(realization%option%mycomm,'test_scatter.out', &
!                              viewer,ierr)
!   call VecView(xx,viewer , ierr)
!   call VecView(field%flow_xx_loc_faces , viewer,ierr)
!   
!   call PetscViewerDestroy(viewer,ierr)
 
 

  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)


  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xz_loc,field%perm_xz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xy_loc,field%perm_xy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yz_loc,field%perm_yz_loc,ONEDOF)


!  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
!         call RichardsUpdateCellPressure(realization)
!  end if  

  
  ! pass #1 for flux terms and accumulation term
   call RichardsResidualPatchMFDLP1(snes,xx,r,realization,ierr)

  ! pass #2 for boundary and source data
   call RichardsResidualPatchMFDLP2(snes,xx,r,realization,ierr)
!      call RichardsCheckMassBalancePatch(realization)

   call VecScatterBegin( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r, &
                                ADD_VALUES,SCATTER_REVERSE, ierr)
   call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces, ierr)

   call VecScatterEnd ( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r,&
                                ADD_VALUES,SCATTER_REVERSE, ierr)


  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RresidualMFD.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RxxMFD.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
 endif


  call PetscLogEventEnd(logging%event_r_residual,ierr)

!  call VecNorm(r, NORM_2, rnorm, ierr)

!  write(*,*) "Residual norm", rnorm

!   write(*,*) "End RichardsResidualMFD"
!   option%test_res = option%test_res + 1
 !  stop   
!   read(*,*) 

#endif
! read(*,*)


end subroutine RichardsResidualMFDLP

! ************************************************************************** !
!
! RichardsResidualPatch1: Computes the interior flux and boundary flux 
!   terms of the residual equation on a single patch
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatch1(snes,xx,r,realization,ierr)

  use water_eos_module

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
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

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
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
  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  call RichardsUpdateAuxVarsPatch(realization)
  patch%aux%Richards%aux_vars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
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
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)

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

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &                  ! distance_gravity = dx*g*n
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(3,iconn))

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      select case (realization%discretization%hydr_flux_method)
        case (TWO_POINT_FLUX)
          call RichardsFlux(rich_aux_vars(ghosted_id_up), &
                            global_aux_vars(ghosted_id_up), &
                            porosity_loc_p(ghosted_id_up), &
                            richards_parameter%sir(1,icap_up), &
                            dd_up,perm_up, &
                            rich_aux_vars(ghosted_id_dn), &
                            global_aux_vars(ghosted_id_dn), &
                            porosity_loc_p(ghosted_id_dn), &
                            richards_parameter%sir(1,icap_dn), &
                            dd_dn,perm_dn, &
                            cur_connection_set%area(iconn), &
                            cur_connection_set%dist(1:3,iconn), &
                            distance_gravity, &
                            upweight,option,v_darcy,Res)
        case (LSM_FLUX)
          call RichardsLSMFlux(rich_aux_vars(ghosted_id_up), &
                               global_aux_vars(ghosted_id_up), &
                               porosity_loc_p(ghosted_id_up), &
                               richards_parameter%sir(1,icap_up), &
                               dd_up,perm_up, &
                               rich_aux_vars(ghosted_id_dn), &
                               global_aux_vars(ghosted_id_dn), &
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
        case default
          option%io_buffer = 'Unknown hydr_flux_method '
          call printErrMsg(option)
      end select

      patch%internal_velocities(1,sum_connection) = v_darcy

#ifdef COMPUTE_INTERNAL_MASS_FLUX
      global_aux_vars(local_id_up)%mass_balance_delta(1,1) = &
        global_aux_vars(local_id_up)%mass_balance_delta(1,1) - Res(1)
#endif

      if (local_id_up>0) then
        r_p(local_id_up) = r_p(local_id_up) + Res(1)
      endif
         
      if (local_id_dn>0) then
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

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*dabs(cur_connection_set%dist(3,iconn))
      icap_dn = patch%sat_func_id(ghosted_id)

      call RichardsBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                rich_aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                richards_parameter%sir(1,icap_dn), &
                                perm_dn, &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(0:3,iconn), &
                                option, &
                                v_darcy,Res)
      patch%boundary_velocities(1,sum_connection) = v_darcy

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_aux_vars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_aux_vars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
        ! contribution to internal 
!        global_aux_vars(ghosted_id)%mass_balance_delta(1) = &
!          global_aux_vars(ghosted_id)%mass_balance_delta(1) + Res(1)
      endif

      r_p(local_id)= r_p(local_id) - Res(1)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)

  !read(*,*) local_id

end subroutine RichardsResidualPatch1

! ************************************************************************** !
!
! RichardsResidualPatch2: Computes the accumulation and source/sink terms of 
!   the residual equation on a single patch
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatch2(snes,xx,r,realization,ierr)

  use water_eos_module

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

  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:)

  PetscReal :: qsrc, qsrc_mol
  PetscReal :: Res(realization%option%nflowdof)


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_ss(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_ss(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_ss => patch%aux%Richards%aux_vars_ss
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

  ! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%steady_state) then
    r_p = r_p - accum_p

    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option,Res) 
      r_p(local_id) = r_p(local_id) + Res(1)
    enddo
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    if(source_sink%flow_condition%rate%itype/=DISTRIBUTED_VOLUMETRIC_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=DISTRIBUTED_MASS_RATE_SS) &
      qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O ! kg/sec -> kmol/sec
        case(SCALED_MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O* & ! kg/sec -> kmol/sec
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_aux_vars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc_mol = qsrc*global_aux_vars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(DISTRIBUTED_VOLUMETRIC_RATE_SS)
          ! qsrc1 = m^3/sec
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* & ! flow = m^3/s
                     global_aux_vars(ghosted_id)%den(1)                  ! den  = kmol/m^3
        case(DISTRIBUTED_MASS_RATE_SS)
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O ! kg/sec -> kmol/sec
      end select
      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
        global_aux_vars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_aux_vars_ss(sum_connection)%mass_balance_delta(1,1) - &
          qsrc_mol
      endif
      r_p(local_id) = r_p(local_id) - qsrc_mol
      ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
      patch%ss_fluid_fluxes(1,sum_connection) = qsrc_mol / &
                                             global_aux_vars(ghosted_id)%den(1)
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%Richards%inactive_cells_exist) then
    do i=1,patch%aux%Richards%n_zero_rows
      r_p(patch%aux%Richards%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

  
#ifdef DASVYAT
      
!   write(*,*) "Internal faces"
!   do i = 1, grid%internal_connection_set_list%first%num_connections
!      write(*,*) "int_flux", i, patch%internal_velocities(option%nphase, i) 
!   end do  

!   write(*,*) "Internal faces"
!   do i = 1, grid%internal_connection_set_list%first%num_connections
!      write(*,*) "int_flux", i, patch%internal_velocities(option%nphase, i) 
!   end do 

!    write(*,*) "Boundary faces"
!   do i = 1, 30
!      write(*,*) "bound_flux", i, patch%boundary_velocities(option%nphase, i) 
!   end do  

!    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_r.out', &
!                              viewer,ierr)   
!    call VecView(realization%field%flow_r,viewer,ierr) 
!    call VecView(r,viewer,ierr) 

!    call PetscViewerDestroy(viewer,ierr)

!     write(*,*) "stop RichardsResidualPatch2"
!     read(*,*)
!     stop

#endif


end subroutine RichardsResidualPatch2



! ************************************************************************** !
!
! RichardsResidualPatchMFD1: Computes flux and accumulation 
!   terms of the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)

  use water_eos_module
  use Logging_module
  use Connection_module
  use Realization_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use MFD_Aux_module
  use MFD_module

  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"     
     

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:),&
                        volume_p(:), xx_loc_faces_p(:), xx_p(:), work_loc_faces_p(:), &
                        perm_xz_loc_p(:), perm_xy_loc_p(:), perm_yz_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:), bc_faces_p(:)


  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn, bound_id 

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  PetscScalar, pointer :: sq_faces(:), e2n_local(:), Smatrix(:,:), face_pr(:)
   PetscScalar, pointer :: neig_den(:), neig_kvr(:), neig_pres(:), bnd(:)
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3), den, ukvr
  PetscInt :: icell, iface, jface, i,j, numfaces, ghost_face_id, ghost_face_jd
  PetscInt :: ghost_neig_id
  PetscViewer :: viewer

  call PetscLogEventBegin(logging%event_flow_residual_mfd1, ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
  discretization => realization%discretization 

!  call RichardsUpdateAuxVarsPatchMFD(realization)
  call RichardsUpdateAuxVarsPatch(realization)


  patch%aux%Richards%aux_vars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call RichardsZeroMassBalDeltaPatch(realization)
  endif


#ifdef DASVYAT_DEBUG
      
!   write(*,*) "Internal faces"
!   do i = 1, grid%internal_connection_set_list%first%num_connections
!      write(*,*) "int_flux", option%myrank, i, patch%internal_velocities(option%nphase, i) 
!   end do  

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_xx_faces_init_mim.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_xx_loc_faces,viewer,ierr) 
!    call VecView(r,viewer,ierr) 

    call PetscViewerDestroy(viewer,ierr)

#endif



! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
!  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr)
 ! call VecGetArrayF90(field%work_loc_faces, work_loc_faces_p, ierr)
 ! call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
 ! call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)
!  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  !print *,' Finished scattering non deriv'


  numfaces = 6 ! hex only
  allocate(sq_faces(numfaces))
  allocate(face_pr(numfaces))
  allocate(neig_den(numfaces))
  allocate(neig_kvr(numfaces))
  allocate(neig_pres(numfaces))
  allocate(bnd(numfaces))
!  allocate(Smatrix(numfaces, numfaces))

  bnd = 0

  !write(*,*) "cell centered"
 ! write(*,*) (xx_p(icell), icell=1,grid%nlmax)  
!  write(*,*) "face centered"
!  write(*,*) (xx_loc_faces_p(iface), iface=1,grid%nlmax_faces)  
!  read(*,*)
!  r_p = 0.

  do icell = 1, grid%nlmax

!   xx_loc_p(icell) = 0

    aux_var => grid%MFD%aux_vars(icell)
    do j = 1, aux_var%numfaces
       ghosted_id = grid%nL2G(icell)
       ghost_face_id = aux_var%face_id_gh(j)
       conn => grid%faces(ghost_face_id)%conn_set_ptr
       jface = grid%faces(ghost_face_id)%id
       sq_faces(j) = conn%area(jface)
       face_pr(j) = xx_loc_faces_p(ghost_face_id)
 
!       if ( (e2n_local(j + (icell-1)*numfaces) == -DIRICHLET_BC).or. &
!                 (e2n_local(j + (icell-1)*numfaces) == -HYDROSTATIC_BC))  then
!             bound_id = grid%fL2B(grid%fG2L(ghost_face_id))
!            write(*,*) "check" ,bound_id, face_pr(j),  bc_faces_p(ghost_face_id)/sq_faces(j), &
!                               face_pr(j) - bc_faces_p(ghost_face_id)/sq_faces(j)
!              face_pr(j) = bc_faces_p(ghost_face_id)/sq_faces(j)
!
!       end if


!       xx_loc_faces_p(ghost_face_id) = conn%cntr(3,jface)     ! dasvyat test
!       xx_loc_p(icell) = xx_loc_p(icell) + conn%cntr(3,jface)/6.
!       write(*,*) conn%cntr(3,jface), xx_loc_faces_p(ghost_face_id)
          if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
               neig_den(j) = global_aux_vars(ghosted_id)%den(1)
          else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
                if (ghosted_id == conn%id_up(jface)) then
                    ghost_neig_id = conn%id_dn(jface)
                else 
                    ghost_neig_id = conn%id_up(jface)
                end if
 
                neig_den(j) = global_aux_vars(ghost_neig_id)%den(1)
          end if
    end do

      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
 

      PermTensor = 0.
      PermTensor(1,1) = perm_xx_loc_p(ghosted_id)
      PermTensor(2,2) = perm_yy_loc_p(ghosted_id)
      PermTensor(3,3) = perm_zz_loc_p(ghosted_id)
      PermTensor(1,3) = perm_xz_loc_p(ghosted_id)
      PermTensor(1,2) = perm_xy_loc_p(ghosted_id)
      PermTensor(2,3) = perm_yz_loc_p(ghosted_id)
      PermTensor(2,1) = PermTensor(1,2) 
      PermTensor(3,1) = PermTensor(1,3) 
      PermTensor(3,2) = PermTensor(2,3) 


        call PetscLogEventBegin(logging%event_flow_flux_mfd, ierr)

      call MFDAuxFluxes(patch, grid, ghosted_id, xx_p(icell:icell), face_pr, aux_var, PermTensor, &
                        rich_aux_vars(ghosted_id), global_aux_vars(ghosted_id), &
                        sq_faces,  bnd, neig_den, neig_kvr, neig_pres, option)
       
        call PetscLogEventEnd(logging%event_flow_flux_mfd, ierr)
      
  end do

!  do iface =1,grid%ngmax_faces
!    write(*,*) iface, "res_p ", r_p(iface), "x", xx_loc_faces_p(iface)
!     write(*,*) xx_loc_faces_p(iface) - work_loc_faces_p(iface)
!  end do

!   write(*,*) "Internal faces", grid%internal_connection_set_list%first%num_connections
!   do iface = 1, grid%internal_connection_set_list%first%num_connections
!     write(*,*) "int_flux",  iface, patch%internal_velocities(option%nphase, iface) 
!   end do  

#ifdef DASVYAT
!   write(*,*) "Boundary faces"
!   do iface = 1, 30 
!      write(*,*) "bound_flux", iface, patch%boundary_velocities(option%nphase, iface)
!   end do  
!   read(*,*)
#endif

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
!  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
!  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
!  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
!  call VecRestoreArrayF90(field%work_loc_faces, work_loc_faces_p, ierr)
!  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
!  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)

  deallocate(sq_faces)
  deallocate(face_pr)
  deallocate(neig_den)
  deallocate(bnd)

#ifdef DASVYAT_DEBUG
  write(*,*) "richards 2822"
  write(*,*) "End RichardsResidualPatchMFD1"
  read(*,*)
!  stop
#endif

  
  call PetscLogEventEnd(logging%event_flow_residual_mfd1, ierr)

#endif

end subroutine RichardsResidualPatchMFD1

! ************************************************************************** !
!
! RichardsResidualPatchMFD2: Computes the boundary and source/sink terms of 
!   the residual equation for MFD discretization
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatchMFD2(snes,xx,r,realization,ierr)

  use water_eos_module
  use Logging_module
  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use MFD_module
  use MFD_aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: i,j
  PetscInt :: local_id, ghosted_id
  


  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:),  accum_p(:), bc_faces_p(:)
  PetscScalar, pointer ::  perm_xx_loc_p(:),  perm_yy_loc_p(:), perm_zz_loc_p(:), flow_xx_p(:)
  PetscScalar, pointer ::  perm_xz_loc_p(:),  perm_xy_loc_p(:), perm_yz_loc_p(:)
  PetscScalar, pointer :: xx_loc_faces_p(:)


  PetscScalar :: qsrc, qsrc_mol


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(coupler_type), pointer :: source_sink, boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection, bc_type, stride


  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
! PetscReal, pointer :: sq_faces(:), rhs(:), bc_g(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: sq_faces(6), rhs(6), bc_g(6), bc_h(6), face_pres(6), bnd(6)
  PetscScalar :: neig_den(6), neig_kvr(6), neig_dkvr_dp(6)
  PetscScalar :: Accum(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscScalar :: Res(realization%option%nflowdof), PermTensor(3,3)
  PetscInt :: icell, iface, jface, numfaces, ghost_face_id, ghost_face_jd, ghost_neig_id
  PetscScalar, pointer :: e2n_local(:)
  


  call PetscLogEventBegin(logging%event_flow_residual_mfd2, ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter

  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_xx, flow_xx_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)


  numfaces = 6 ! hex only
#if 0  
  allocate(sq_faces(numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"sq_faces allocation error" 

  allocate(face_pres(numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"sq_faces allocation error" 

  allocate(rhs(numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"rhs allocation error" 

  allocate(bc_g(numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"bc_g allocation error" 

  allocate(bc_h(numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"bc_h allocation error" 

  allocate(bnd (numfaces), stat = ierr)
  if (ierr /= 0) write(*,*)"bnd allocation error" 
#endif
  stride = 6 !hex only

  ! open(unit=36, FILE = "rhs.dat", STATUS = "REPLACE")

  do local_id = 1, grid%nlmax

        ghosted_id = grid%nL2G(local_id)
        bc_g = 0
        bc_h = 0
        bnd = 0
 

        aux_var => grid%MFD%aux_vars(local_id)
        do j = 1, aux_var%numfaces
          ghost_face_id = aux_var%face_id_gh(j)
          conn => grid%faces(ghost_face_id)%conn_set_ptr
          jface = grid%faces(ghost_face_id)%id
          sq_faces(j) = conn%area(jface)
          face_pres(j) = xx_loc_faces_p(ghost_face_id)
         
  !       write(*,*) "Face ", conn%itype, jface
 
  !        write(*,*) j, e2n_local(j + (local_id-1)*stride)
          if ( (e2n_local(j + (local_id-1)*stride) == -DIRICHLET_BC).or. &
                 (e2n_local(j + (local_id-1)*stride) == -HYDROSTATIC_BC))  then
            bc_g(j) = bc_faces_p(ghost_face_id)
            face_pres(j) = 0.
            bnd(j) = 1
          else if ( e2n_local(j + (local_id-1)*stride) == -NEUMANN_BC ) then
 !           write(*,*) "Neumann", bc_faces_p(ghost_face_id)
            bc_h(j) = bc_faces_p(ghost_face_id)
          else if (e2n_local(j + (local_id-1)*stride) < 0) then 
            call printMsg(option,'Type of BC is not supported by MFD mode')
            stop
          end if

          if (conn%itype == BOUNDARY_CONNECTION_TYPE) then

!               neig_den(j) = global_aux_vars(ghosted_id)%den(1)
               neig_den(j) = global_aux_vars_bc(jface)%den(1) 

#ifdef USE_ANISOTROPIC_MOBILITY
               neig_kvr(j) = rich_aux_vars(ghosted_id)%kvr_x
               neig_dkvr_dp(j) = rich_aux_vars(ghosted_id)%dkvr_x_dp
#else
!               neig_kvr(j) = rich_aux_vars(ghosted_id)%kvr
               neig_kvr(j) = rich_aux_vars_bc(jface)%kvr
!               neig_dkvr_dp(j) = rich_aux_vars(ghosted_id)%dkvr_dp
               neig_dkvr_dp(j) = 0.
#endif

          else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
                if (ghosted_id == conn%id_up(jface)) then
                    ghost_neig_id = conn%id_dn(jface)
                else 
                    ghost_neig_id = conn%id_up(jface)
                end if
 
                neig_den(j) = global_aux_vars(ghost_neig_id)%den(1)
#ifdef USE_ANISOTROPIC_MOBILITY
                neig_kvr(j) = rich_aux_vars(ghost_neig_id)%kvr_x 
                neig_dkvr_dp(j) = rich_aux_vars(ghost_neig_id)%dkvr_x_dp
#else
                neig_kvr(j) = rich_aux_vars(ghost_neig_id)%kvr
                neig_dkvr_dp(j) = rich_aux_vars(ghost_neig_id)%dkvr_dp
#endif
          end if
 

        end do
!          write(*,*) "neig_den", (neig_den(i),i=1,6)
         

        !geh - Ignore inactive cells with inactive materials
        if (patch%imat(ghosted_id) <= 0) cycle
        call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option, Res)

            

        Accum(1) = Res(1) - accum_p(local_id)
!        write(*,*) "Accum ", Accum(1)
!		Accum(1) = 0
!        write(*,*) "accum ",accum_p(local_id), "Res ", Res(1),  "diff ", Accum(1)
!       write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)
 
        source_f = 0.
        PermTensor = 0.
        PermTensor(1,1) = perm_xx_loc_p(ghosted_id)
        PermTensor(2,2) = perm_yy_loc_p(ghosted_id)
        PermTensor(3,3) = perm_zz_loc_p(ghosted_id)
        PermTensor(1,3) = perm_xz_loc_p(ghosted_id)
        PermTensor(1,2) = perm_xy_loc_p(ghosted_id)
        PermTensor(2,3) = perm_yz_loc_p(ghosted_id)
        PermTensor(3,1) = PermTensor(1,3)
        PermTensor(2,1) = PermTensor(1,2)
        PermTensor(3,2) = PermTensor(2,3)

        call PetscLogEventBegin(logging%event_flow_rhs_mfd, ierr)

         

!        call MFDAuxGenerateRhs(patch, grid, ghosted_id, PermTensor, bc_g, source_f, bc_h, aux_var, &
!                                 rich_aux_vars(ghosted_id),&
!                                 global_aux_vars(ghosted_id),&
!                                 Accum, &
!                                 porosity_loc_p(ghosted_id), volume_p(local_id),&
!                                 flow_xx_p(local_id:local_id), face_pres, bnd,&                                 
!                                 sq_faces, neig_den, neig_kvr, neig_dkvr_dp, option, rhs) 

        call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr)
 
        call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr)

   !     stop

   !     stop
     
        do iface=1, aux_var%numfaces

          ghost_face_id = aux_var%face_id_gh(iface) 

          if (e2n_local((local_id -1)*numfaces + iface) < 0) then
            if ((e2n_local((local_id -1)*numfaces + iface) == -DIRICHLET_BC).or. &
                (e2n_local((local_id -1)*numfaces + iface) == -HYDROSTATIC_BC).or.&
                (e2n_local((local_id -1)*numfaces + iface) == -SEEPAGE_BC).or.&
                (e2n_local((local_id -1)*numfaces + iface) == -CONDUCTANCE_BC)) then

                r_p(ghost_face_id) = 0.
            else 
               r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface) 
            end if
          else
                r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface)
          end if
        end do
    
!          write(*,*) (rhs(iface),iface=1,6)
!          write(*,*)
       
  enddo


!  close(36)

!  open(unit=35, FILE = "res.dat", STATUS = "REPLACE")
!  do iface =1,grid%ngmax_faces
!    write(35,*) "residual_p ", grid%fG2P(iface), r_p(iface), xx_loc_faces_p(iface)
!  end do
!   close(35)

!  write(*,*)
!  read(*,*)

#if 0 
  deallocate(sq_faces)
  deallocate(rhs)
  deallocate(bc_g)
  deallocate(bc_h)
  deallocate(face_pres)
  deallocate(bnd)
#endif

  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx, flow_xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)

!  write(*,*) "End RichardsResidualPatchMFD2"
!  read(*,*)

!  stop 
  call PetscLogEventEnd(logging%event_flow_residual_mfd2, ierr)

#endif


end subroutine RichardsResidualPatchMFD2


! ************************************************************************** !
!
! RichardsResidualPatchMFD1: Computes flux and accumulation 
!   terms of the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatchMFDLP1(snes,xx,r,realization,ierr)

  use water_eos_module
  use Logging_module
  use Connection_module
  use Realization_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use MFD_Aux_module
  use MFD_module

  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"     
     

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:),&
                        volume_p(:), xx_loc_faces_p(:), xx_p(:), work_loc_faces_p(:), &
                        perm_xz_loc_p(:), perm_xy_loc_p(:), perm_yz_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:), bc_faces_p(:)


  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:) 
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(richards_auxvar_type) :: test_rich_aux_vars
  type(global_auxvar_type) :: test_global_aux_vars
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn, bound_id, LP_cell_id 

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  PetscScalar, pointer :: sq_faces(:), e2n_local(:), Smatrix(:,:), face_pr(:)
  PetscScalar, pointer :: neig_den(:), neig_ukvr(:), neig_pres(:), bnd(:), bc_h(:)
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3), den, ukvr
  PetscInt :: icell, iface, jface, i,j, numfaces, ghost_face_id, ghost_face_jd
  PetscInt :: ghost_neig_id
  PetscViewer :: viewer

  call PetscLogEventBegin(logging%event_flow_residual_mfd1, ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  discretization => realization%discretization 

  call RichardsUpdateAuxVarsPatchMFDLP(realization)


  patch%aux%Richards%aux_vars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call RichardsZeroMassBalDeltaPatch(realization)
  endif

!  write(*,*) "RichardsResidualPatchMFDLP1"
!  read(*,*) 


! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)


  numfaces = 6 ! hex only
  allocate(sq_faces(numfaces))
  allocate(face_pr(numfaces))
  allocate(neig_den(numfaces))
  allocate(neig_ukvr(numfaces))
  allocate(neig_pres(numfaces))
  allocate(bnd(numfaces))
  allocate(bc_h(numfaces))



  !write(*,*) "cell centered"
 ! write(*,*) (xx_p(icell), icell=1,grid%nlmax)  
!  write(*,*) "face centered"
!  write(*,*) (xx_loc_faces_p(iface), iface=1,grid%nlmax_faces)  
!  read(*,*)
!  r_p = 0.

  do icell = 1, grid%nlmax

    bnd = 0;
!   xx_loc_p(icell) = 0

    aux_var => grid%MFD%aux_vars(icell)


    do j = 1, aux_var%numfaces
       ghosted_id = grid%nL2G(icell)
       ghost_face_id = aux_var%face_id_gh(j)
       conn => grid%faces(ghost_face_id)%conn_set_ptr
       jface = grid%faces(ghost_face_id)%id
       sq_faces(j) = conn%area(jface)
       face_pr(j) = xx_loc_faces_p(ghost_face_id)
 
       if ( (e2n_local(j + (icell-1)*numfaces) == -1.0*DIRICHLET_BC).or. &
                (e2n_local(j + (icell-1)*numfaces) == -1.0*HYDROSTATIC_BC))  then
           bnd(j) = 1
       else if ( e2n_local(j + (icell-1)*numfaces) == -1.0*NEUMANN_BC ) then
 !          write(*,*) "Neumann", bc_faces_p(ghost_face_id)
            bnd(j) = 2
       else if (e2n_local(j + (icell-1)*numfaces)< 0) then 
           call printMsg(option,'Type of BC is not supported by MFD mode')
           stop
       end if


!       xx_loc_faces_p(ghost_face_id) = conn%cntr(3,jface)     ! dasvyat test
!       xx_loc_p(icell) = xx_loc_p(icell) + conn%cntr(3,jface)/6.
!       write(*,*) conn%cntr(3,jface), xx_loc_faces_p(ghost_face_id)
          if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
                  neig_pres(j) = xx_loc_faces_p(ghost_face_id)

               if (bnd(j) > 0) then
                  call RichardsAuxVarInit(test_rich_aux_vars, option)  
                  call GlobalAuxVarInit(test_global_aux_vars, option)
!
                  call RichardsAuxVarCompute(neig_pres(j:j), test_rich_aux_vars, &
                       test_global_aux_vars, &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                       option)
               end if
               if (bnd(j)==1) then


                  neig_den(j) = test_global_aux_vars%den(1)
!                  write(*,*) "jface ", jface, neig_pres(j)

!                  neig_den(j) = global_aux_vars_bc(jface)%den(1)


!                  neig_kvr(j) = rich_aux_vars_bc(jface)%kvr
                  neig_ukvr(j) = test_rich_aux_vars%kvr
!                  neig_dkvr_dp(j) =  rich_aux_vars(ghosted_id)%dkvr_dp
!                  neig_dkvr_dp(j) =  test_rich_aux_vars%dkvr_dp
               else if (bnd(j)==2) then
!                  if (bc_faces_p(ghost_face_id) > 0) then
!                     neig_den(j) = test_global_aux_vars%den(1)
!                  else
                     neig_den(j) = global_aux_vars(ghosted_id)%den(1)
                     
 !                 end if
                else
                  neig_den(j) = global_aux_vars(ghosted_id)%den(1)
                  neig_ukvr(j) = rich_aux_vars(ghosted_id)%kvr
                end if
                  if (bnd(j) > 0) call GlobalAuxVarStrip(test_global_aux_vars) 
 

          else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
                if (ghosted_id == conn%id_up(jface)) then
                    ghost_neig_id = conn%id_dn(jface)
                else 
                    ghost_neig_id = conn%id_up(jface)
                end if
 
                neig_den(j) = global_aux_vars(ghost_neig_id)%den(1)
                neig_pres(j) = xx_loc_faces_p(grid%ngmax_faces + ghost_neig_id)
#ifdef USE_ANISOTROPIC_MOBILITY
                neig_ukvr(j) = rich_aux_vars(ghost_neig_id)%kvr_x
#else
                neig_ukvr(j) = rich_aux_vars(ghost_neig_id)%kvr
#endif
          end if
    end do

      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
 

      PermTensor = 0.
      PermTensor(1,1) = perm_xx_loc_p(ghosted_id)
      PermTensor(2,2) = perm_yy_loc_p(ghosted_id)
      PermTensor(3,3) = perm_zz_loc_p(ghosted_id)
      PermTensor(1,3) = perm_xz_loc_p(ghosted_id)
      PermTensor(1,2) = perm_xy_loc_p(ghosted_id)
      PermTensor(2,3) = perm_yz_loc_p(ghosted_id)
      PermTensor(2,1) = PermTensor(1,2) 
      PermTensor(3,1) = PermTensor(1,3) 
      PermTensor(3,2) = PermTensor(2,3) 


      call PetscLogEventBegin(logging%event_flow_flux_mfd, ierr)

      LP_cell_id = grid%ngmax_faces + ghosted_id

      call MFDAuxFluxes(patch, grid, ghosted_id, xx_loc_faces_p(LP_cell_id:LP_cell_id), face_pr, aux_var, PermTensor, &
                        rich_aux_vars(ghosted_id), global_aux_vars(ghosted_id), &
                        sq_faces,  bnd, neig_den, neig_ukvr, neig_pres, option)
       
      call PetscLogEventEnd(logging%event_flow_flux_mfd, ierr)
      
  end do

!  do iface =1,grid%ngmax_faces
!    write(*,*) iface, "res_p ", r_p(iface), "x", xx_loc_faces_p(iface)
!     write(*,*) xx_loc_faces_p(iface) - work_loc_faces_p(iface)
!  end do

!   write(*,*) "Internal faces", grid%internal_connection_set_list%first%num_connections
!   do iface = 1, grid%internal_connection_set_list%first%num_connections
!     write(*,*) "int_flux",  iface, patch%internal_velocities(option%nphase, iface) 
!   end do  

#ifdef DASVYAT
!   write(*,*) "Boundary faces"
!   do iface = 1, 30 
!      write(*,*) "bound_flux", iface, patch%boundary_velocities(option%nphase, iface)
!   end do  
!   read(*,*)
#endif

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)

  deallocate(sq_faces)
  deallocate(face_pr)
  deallocate(neig_den)
  deallocate(neig_ukvr)
  deallocate(neig_pres)
  deallocate(bnd)

#ifdef DASVYAT_DEBUG
  write(*,*) "richards 2822"
  write(*,*) "End RichardsResidualPatchMFD1"
!  read(*,*)
!  stop
#endif

  
  call PetscLogEventEnd(logging%event_flow_residual_mfd1, ierr)

#endif

end subroutine RichardsResidualPatchMFDLP1

! ************************************************************************** !
!
! RichardsResidualPatchMFDLP2: Computes the boundary and source/sink terms of 
!   the residual equation for MFD discretization
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatchMFDLP2(snes,xx,r,realization,ierr)

  use water_eos_module
  use Logging_module
  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use MFD_module
  use MFD_aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: i,j
  PetscInt :: local_id, ghosted_id
  


  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:),  accum_p(:), bc_faces_p(:)
  PetscScalar, pointer ::  perm_xx_loc_p(:),  perm_yy_loc_p(:), perm_zz_loc_p(:), flow_xx_p(:)
  PetscScalar, pointer ::  perm_xz_loc_p(:),  perm_xy_loc_p(:), perm_yz_loc_p(:)
  PetscScalar, pointer :: xx_loc_faces_p(:)


  PetscScalar :: qsrc, qsrc_mol, tempreal


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(richards_auxvar_type) :: test_rich_aux_vars
  type(global_auxvar_type) :: test_global_aux_vars
  type(coupler_type), pointer :: source_sink, boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection, bc_type, stride


  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
! PetscReal, pointer :: sq_faces(:), rhs(:), bc_g(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: sq_faces(6), rhs(7), bc_g(6), bc_h(6), face_pres(6), bnd(6)
  PetscScalar :: neig_den(6), neig_kvr(6), neig_dkvr_dp(6), neig_pres(6)
  PetscScalar :: Accum(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscScalar :: Res(realization%option%nflowdof), PermTensor(3,3)
  PetscInt :: icell, iface, jface, numfaces, ghost_face_id, ghost_face_jd, ghost_neig_id, LP_cell_id
  PetscScalar, pointer :: e2n_local(:)
  


  call PetscLogEventBegin(logging%event_flow_residual_mfd2, ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter

  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc



! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_xx, flow_xx_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)



  numfaces = 6 ! hex only
  stride = 6 !hex only

  ! open(unit=36, FILE = "rhs.dat", STATUS = "REPLACE")

  do local_id = 1, grid%nlmax

        ghosted_id = grid%nL2G(local_id)
        bc_g = 0
        bc_h = 0
        bnd = 0
 

        aux_var => grid%MFD%aux_vars(local_id)
        do j = 1, aux_var%numfaces
          ghost_face_id = aux_var%face_id_gh(j)
          conn => grid%faces(ghost_face_id)%conn_set_ptr
          jface = grid%faces(ghost_face_id)%id
          sq_faces(j) = conn%area(jface)
          face_pres(j) = xx_loc_faces_p(ghost_face_id)
         
 
          if ( (e2n_local(j + (local_id-1)*stride) == -DIRICHLET_BC).or. &
                 (e2n_local(j + (local_id-1)*stride) == -HYDROSTATIC_BC))  then
            bc_g(j) = bc_faces_p(ghost_face_id)
            face_pres(j) = 0.
            bnd(j) = 1
          else if ( e2n_local(j + (local_id-1)*stride) == -NEUMANN_BC ) then
 !           write(*,*) "Neumann", bc_faces_p(ghost_face_id)
            bc_h(j) = bc_faces_p(ghost_face_id)
            bnd(j) = 2
          else if (e2n_local(j + (local_id-1)*stride) < 0) then 
            call printMsg(option,'Type of BC is not supported by MFD mode')
            stop
          end if

          if (conn%itype == BOUNDARY_CONNECTION_TYPE) then

               neig_pres(j) = xx_loc_faces_p(ghost_face_id)

               if (bnd(j) > 0) then
                  call RichardsAuxVarInit(test_rich_aux_vars, option)  
                  call GlobalAuxVarInit(test_global_aux_vars, option)
!
                  call RichardsAuxVarCompute(neig_pres(j:j), test_rich_aux_vars, &
                       test_global_aux_vars, &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                       option)
               end if
               if (bnd(j)==1) then


                  neig_den(j) = test_global_aux_vars%den(1)
!                  write(*,*) "jface ", jface, neig_pres(j)

!                  neig_den(j) = global_aux_vars_bc(jface)%den(1)


!                  neig_kvr(j) = rich_aux_vars_bc(jface)%kvr
                  neig_kvr(j) = test_rich_aux_vars%kvr
!                  neig_dkvr_dp(j) =  rich_aux_vars(ghosted_id)%dkvr_dp
!                  neig_dkvr_dp(j) =  test_rich_aux_vars%dkvr_dp
                  neig_dkvr_dp(j) =  0.
               else if (bnd(j)==2) then
!                  if (bc_h(j) > 0) then
!                     neig_den(j) = test_global_aux_vars%den(1)
!                  else
                     neig_den(j) = global_aux_vars(ghosted_id)%den(1)
                     
!                  end if
                else
                  neig_den(j) = global_aux_vars(ghosted_id)%den(1)
                  neig_dkvr_dp(j) = rich_aux_vars(ghosted_id)%dkvr_dp
                  neig_kvr(j) = rich_aux_vars(ghosted_id)%kvr
                end if
                  if (bnd(j) > 0) call GlobalAuxVarStrip(test_global_aux_vars) 

          else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
                if (ghosted_id == conn%id_up(jface)) then
                    ghost_neig_id = conn%id_dn(jface)
                else 
                    ghost_neig_id = conn%id_up(jface)
                end if
 
                neig_den(j) = global_aux_vars(ghost_neig_id)%den(1)
                neig_kvr(j) = rich_aux_vars(ghost_neig_id)%kvr 
                neig_dkvr_dp(j) = rich_aux_vars(ghost_neig_id)%dkvr_dp
                neig_pres(j) = xx_loc_faces_p(grid%ngmax_faces + ghost_neig_id)
          end if
 

        end do
         

        !geh - Ignore inactive cells with inactive materials
        if (patch%imat(ghosted_id) <= 0) cycle
        call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option, Res)

            

        Accum(1) = Res(1) - accum_p(local_id)
!        write(*,*) "Accum ", Accum(1)
!		Accum(1) = 0
!        if (ghosted_id==43) then
!            write(*,*) "accum ",accum_p(local_id), "Res ", Res(1),  "diff ", Accum(1)
!            write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)
!        end if
 
        source_f = 0.
        PermTensor = 0.
        PermTensor(1,1) = perm_xx_loc_p(ghosted_id)
        PermTensor(2,2) = perm_yy_loc_p(ghosted_id)
        PermTensor(3,3) = perm_zz_loc_p(ghosted_id)
        PermTensor(1,3) = perm_xz_loc_p(ghosted_id)
        PermTensor(1,2) = perm_xy_loc_p(ghosted_id)
        PermTensor(2,3) = perm_yz_loc_p(ghosted_id)
        PermTensor(3,1) = PermTensor(1,3)
        PermTensor(2,1) = PermTensor(1,2)
        PermTensor(3,2) = PermTensor(2,3)

        call PetscLogEventBegin(logging%event_flow_rhs_mfd, ierr)

         

        LP_cell_id = grid%ngmax_faces + ghosted_id
         
!        if (option%myrank == 0) then  
!              write(*,*) LP_cell_id,  xx_loc_faces_p(LP_cell_id)
!        end if  

        call MFDAuxGenerateRhs_LP(patch, grid, ghosted_id, PermTensor, bc_g, source_f, bc_h, aux_var, &
                                 rich_aux_vars(ghosted_id),&
                                 global_aux_vars(ghosted_id),&
                                 Accum, &
                                 porosity_loc_p(ghosted_id), volume_p(local_id),&
                                 xx_loc_faces_p(LP_cell_id:LP_cell_id), face_pres, bnd,&                                 
                                 sq_faces, neig_den, neig_kvr, neig_dkvr_dp, neig_pres, option, rhs) 


        call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr)
 

     
        do iface=1, aux_var%numfaces

          ghost_face_id = aux_var%face_id_gh(iface) 
          conn => grid%faces(ghost_face_id)%conn_set_ptr
          jface = grid%faces(ghost_face_id)%id

          if (e2n_local((local_id -1)*numfaces + iface) < 0) then
            if ((e2n_local((local_id -1)*numfaces + iface) == -DIRICHLET_BC).or. &
                (e2n_local((local_id -1)*numfaces + iface) == -HYDROSTATIC_BC).or.&
                (e2n_local((local_id -1)*numfaces + iface) == -SEEPAGE_BC).or.&
                (e2n_local((local_id -1)*numfaces + iface) == -CONDUCTANCE_BC)) then

                r_p(ghost_face_id) = 0.
            else 
               r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface) 
            end if
          else
                r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface)
          end if
        end do
   
 !       if (grid%
 
        r_p(grid%ngmax_faces + ghosted_id) = r_p(grid%ngmax_faces + ghosted_id) + rhs(numfaces + 1)

!        if (option%myrank == 0) then  
!       if (option%myrank== 0.and.ghosted_id==10) then
!        if (grid%x(ghosted_id)<=300.and.grid%y(ghosted_id)<100) then
!           write(*,*) "Coordinates:", grid%x(ghosted_id), grid%y(ghosted_id), grid%z(ghosted_id)
!           write(*,*) xx_loc_faces_p(LP_cell_id)
!           write(*,*) option%myrank, (face_pres(iface),iface=1,6)
!        end if
  enddo

!  close(36)

!    open(unit=35, FILE = "res.dat", STATUS = "REPLACE")
!    do iface =1,grid%ngmax_faces + grid%ngmax
!      if (option%myrank == 0) then
!         write(35,*) "residual_p ",  iface, r_p(iface), xx_loc_faces_p(iface)
!      else if (option%myrank == 1) then
!         write(36,*) "residual_p ",  iface, r_p(iface), xx_loc_faces_p(iface)
!      end if
!    end do
!    close(35)

!  write(*,*)
!  read(*,*)

  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx, flow_xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)



!  write(*,*) "End RichardsResidualPatchMFDLP2"
!  read(*,*)

!   stop 
  call PetscLogEventEnd(logging%event_flow_residual_mfd2, ierr)

#endif


end subroutine RichardsResidualPatchMFDLP2



! ************************************************************************** !
!
! RichardsJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsJacobian(snes,xx,A,B,flag,realization,ierr)

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
!
! RichardsJacobian: Computes the Jacobian for MFD Discretization
! author: Daniil Svyatskiy
! date: 09/17/10
!
! ************************************************************************** !
subroutine RichardsJacobianMFD(snes,xx,A,B,flag,realization,ierr)

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

  call RichardsJacobianPatchMFD(snes,xx,J,J,flag,realization,ierr)

  if (realization%debug%matview_Jacobian) then
#if 0  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
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

  call PetscLogEventEnd(logging%event_r_jacobian,ierr)
  
end subroutine RichardsJacobianMFD


subroutine RichardsJacobianMFDLP(snes,xx,A,B,flag,realization,ierr)

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

  call RichardsJacobianPatchMFDLP(snes,xx,J,J,flag,realization,ierr)

  if (realization%debug%matview_Jacobian) then
#if 1  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
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

  call PetscLogEventEnd(logging%event_r_jacobian,ierr)

#if 0
    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx_faces.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_dxx_faces,viewer,ierr) 

    call PetscViewerDestroy(viewer,ierr)
 

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy_faces.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_yy_faces,viewer,ierr) 
    call PetscViewerDestroy(viewer,ierr)
#endif


!  write(*,*) "Exit RichardsJacobianMFDLP"
!  read(*,*)
!  stop

end subroutine RichardsJacobianMFDLP
                
! ************************************************************************** !
!
! RichardsJacobianPatch1: Computes the interior flux and boundary flux 
!   terms of the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsJacobianPatch1(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

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

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscInt :: icap_up,icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:) 
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  PetscInt, pointer :: cell_neighbors(:,:)
  
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

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

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  
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
      perm_up = perm_xx_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*dabs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*dabs(cur_connection_set%dist(3,iconn))
    
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      select case (realization%discretization%hydr_flux_method)
        case (TWO_POINT_FLUX)
          call RichardsFluxDerivative(rich_aux_vars(ghosted_id_up), &
                                      global_aux_vars(ghosted_id_up), &
                                      porosity_loc_p(ghosted_id_up), &
                                      richards_parameter%sir(1,icap_up), &
                                      dd_up,perm_up, &
                                      rich_aux_vars(ghosted_id_dn), &
                                      global_aux_vars(ghosted_id_dn), &
                                      porosity_loc_p(ghosted_id_dn), &
                                      richards_parameter%sir(1,icap_dn), &
                                      dd_dn,perm_dn, &
                                      cur_connection_set%area(iconn), &
                                      cur_connection_set%dist(1:3,iconn),&
                                      distance_gravity, &
                                      upweight,option,&
                                      patch%saturation_function_array(icap_up)%ptr,&
                                      patch%saturation_function_array(icap_dn)%ptr,&
                                      Jup,Jdn)
        case (LSM_FLUX)
          call RichardsLSMFluxDerivative(rich_aux_vars(ghosted_id_up), &
                                         global_aux_vars(ghosted_id_up), &
                                         porosity_loc_p(ghosted_id_up), &
                                         richards_parameter%sir(1,icap_up), &
                                         dd_up,perm_up, &
                                         rich_aux_vars(ghosted_id_dn), &
                                         global_aux_vars(ghosted_id_dn), &
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
        case default
          option%io_buffer = 'Unknown hydr_flux_method '
          call printErrMsg(option)
      end select

      if (local_id_up > 0) then
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

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*dabs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*dabs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*dabs(cur_connection_set%dist(3,iconn))
      icap_dn = patch%sat_func_id(ghosted_id) 

      call RichardsBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                rich_aux_vars_bc(sum_connection), &
                                global_aux_vars_bc(sum_connection), &
                                rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                richards_parameter%sir(1,icap_dn), &
                                perm_dn, &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(0:3,iconn), &
                                option, &
                                patch%saturation_function_array(icap_dn)%ptr,&
                                Jdn)
      Jdn = -Jdn
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
  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)

end subroutine RichardsJacobianPatch1

! ************************************************************************** !
!
! RichardsJacobianPatch2: Computes the accumulation and source/sink terms of 
!   the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsJacobianPatch2(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

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

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:)
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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  PetscInt :: flow_pc
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  
  if (.not.option%steady_state) then
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    icap = patch%sat_func_id(ghosted_id)
    call RichardsAccumDerivative(rich_aux_vars(ghosted_id), &
                              global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option, &
                              patch%saturation_function_array(icap)%ptr,&
                              Jup) 
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
    
    if(source_sink%flow_condition%rate%itype/=DISTRIBUTED_VOLUMETRIC_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=DISTRIBUTED_MASS_RATE_SS) &
      qsrc = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS,SCALED_MASS_RATE_SS,DISTRIBUTED_MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_aux_vars(ghosted_id)%dden_dp*FMWH2O
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_aux_vars(ghosted_id)%dden_dp*FMWH2O* &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(DISTRIBUTED_VOLUMETRIC_RATE_SS)
          Jup(1,1) = -source_sink%flow_aux_real_var(ONE_INTEGER,iconn)* &
                    rich_aux_vars(ghosted_id)%dden_dp*FMWH2O

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
  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

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
!
! RichardsJacobianPatch1: Computes local condensed matrices
!   for the Jacobian
! author: Daniil Svyatskiy
! date: 09/17/10
!
! ************************************************************************** !
subroutine RichardsJacobianPatchMFD (snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module
  use mfd_aux_module
  use Connection_module
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use MFD_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr


#ifdef DASVYAT

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
                        perm_xz_loc_p(:), perm_xy_loc_p(:), perm_yz_loc_p(:)
  PetscReal, pointer :: accum_p(:), xx_p(:)
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal, pointer :: J(:)       
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, iface, jface, icell, numfaces
  PetscInt :: sum_connection  
  PetscReal :: ukvr, den, dden_dp, dkr_dp, dp_dlmd, dkvr_x_dp, PermTensor(3,3) 
  PetscInt, pointer :: ghosted_face_id(:), bound_id(:)
  PetscReal, pointer :: bc_faces_p(:), e2n_local(:)
  PetscReal, pointer :: face_pres(:), xx_faces_p(:), volume_p(:), sq_faces(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars


!  write(*,*) "ENTER MFD JACOBIAN"


  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_xx_loc, xx_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)

     numfaces = 6

     allocate(J(numfaces*numfaces))
     allocate(ghosted_face_id(numfaces))
     allocate(face_pres(numfaces))
     allocate(sq_faces(numfaces))
     allocate(bound_id(numfaces))


  do local_id = 1, grid%nlmax
     aux_var => grid%MFD%aux_vars(local_id)
!     numfaces = aux_var%numfaces

     ghosted_id = grid%nL2G(local_id)
#ifdef USE_ANISOTROPIC_MOBILITY         
     ukvr = rich_aux_vars(ghosted_id)%kvr_x 
     dkvr_x_dp = rich_aux_vars(ghosted_id)%dkvr_x_dp
#else
     ukvr = rich_aux_vars(ghosted_id)%kvr
     dkvr_x_dp = rich_aux_vars(ghosted_id)%dkvr_dp
#endif
     dden_dp = rich_aux_vars(ghosted_id)%dden_dp
     den =  global_aux_vars(ghosted_id)%den(1)

     bound_id = 0

     do iface = 1, numfaces
        ghosted_face_id(iface) = aux_var%face_id_gh(iface) - 1
        face_pres(iface) = xx_faces_p(aux_var%face_id_gh(iface))

        conn => grid%faces(ghosted_face_id(iface)+1)%conn_set_ptr
        jface = grid%faces(ghosted_face_id(iface)+1)%id
        sq_faces(iface) = conn%area(jface)

        

        if ( (e2n_local(iface + (local_id-1)*numfaces) == -DIRICHLET_BC).or. &
                 (e2n_local(iface + (local_id-1)*numfaces) == -HYDROSTATIC_BC))  then
            bound_id(iface) = 1
        end if

     end do


   PermTensor = 0.
   PermTensor(1,1) = perm_xx_loc_p(ghosted_id) 
   PermTensor(2,2) = perm_yy_loc_p(ghosted_id) 
   PermTensor(3,3) = perm_zz_loc_p(ghosted_id) 
   PermTensor(1,3) = perm_xz_loc_p(ghosted_id) 
   PermTensor(1,2) = perm_xy_loc_p(ghosted_id) 
   PermTensor(2,3) = perm_yz_loc_p(ghosted_id) 
   PermTensor(3,1) = PermTensor(1,3)    
   PermTensor(2,1) = PermTensor(1,2)    
   PermTensor(3,2) = PermTensor(2,3)    

    J = 0.

   call MFDAuxJacobianLocal( grid, aux_var, &
                                       rich_aux_vars(ghosted_id), global_aux_vars(ghosted_id), &
                                       sq_faces, option, J)

   do iface = 1, numfaces
     if (bound_id(iface) == 1) then
       do jface = 1, numfaces
          J((iface-1)*numfaces + jface) = 0.
          J((jface-1)*numfaces + iface) = 0.
       end do
       J((iface-1)*numfaces + iface) = 1.
     end if
   end do


!     write(*,*) "Stiff Matrix"

!     do iface = 1, numfaces
!          write(*,*) (J((iface-1)*numfaces + jface), jface=1,numfaces)
!     end do


!    read(*,*)


     call MatSetValuesLocal(A, numfaces, ghosted_face_id, numfaces, ghosted_face_id, &
                                        J, ADD_VALUES,ierr)
!     call MatSetValuesLocal(A, 1, ghosted_face_id(1), 1, ghosted_face_id(1), &
!                                        1d0, ADD_VALUES,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do



     deallocate(ghosted_face_id)
     deallocate(J)
     deallocate(face_pres)
     deallocate(sq_faces)
     deallocate(bound_id)

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
        

!          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
!                                         Jup,ADD_VALUES,ierr)
!          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &


  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_mfd.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xz_loc, perm_xz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xy_loc, perm_xy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yz_loc, perm_yz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)

#endif

#ifdef DASVYAT
!  write(*,*) "EXIT MFD JACOBIAN"
!  write(*,*) "richard 4587"
!  stop
#endif

end subroutine RichardsJacobianPatchMFD

! ************************************************************************** !
!
! RichardsJacobianPatch1: Computes local condensed matrices
!   for the Jacobian
! author: Daniil Svyatskiy
! date: 09/17/10
!
! ************************************************************************** !
subroutine RichardsJacobianPatchMFDLP (snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module
  use mfd_aux_module
  use Connection_module
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use MFD_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr


#ifdef DASVYAT

  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal, pointer :: J(:), Jbl(:)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, iface, jface, icell, numfaces, cell_LP_id
  PetscInt :: sum_connection, nrow, ncol  
  PetscInt, pointer :: bound_id(:), ghosted_LP_id(:), neig_LP_id(:)
  PetscReal, pointer :: e2n_local(:),  sq_faces(:), test_p(:)
  PetscReal :: diag
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars
 

  call VecGetArrayF90(grid%e2n, e2n_local, ierr)

     numfaces = 6

    allocate(J((numfaces+1)*(numfaces+1)))
!    allocate(Jbl((2*numfaces+1)*(numfaces+1)))
    allocate(bound_id(numfaces))
    allocate(sq_faces(numfaces))
    allocate(ghosted_LP_id(numfaces+1+numfaces))
    allocate(neig_LP_id(numfaces))

    ncol = 2*numfaces+1
    nrow = numfaces+1

  do local_id = 1, grid%nlmax
     aux_var => grid%MFD%aux_vars(local_id)
!     numfaces = aux_var%numfaces
!     write(*,*) "RichardsJacobianPatchMFDLP", local_id
     ghosted_id = grid%nL2G(local_id)
     cell_LP_id = grid%ngmax_faces + ghosted_id - 1

     bound_id = 0

    do iface = 1, numfaces
        ghosted_LP_id(iface) = aux_var%face_id_gh(iface) - 1

        conn => grid%faces(ghosted_LP_id(iface)+1)%conn_set_ptr
        jface = grid%faces(ghosted_LP_id(iface)+1)%id
        sq_faces(iface) = conn%area(jface)

         if ( (e2n_local(iface + (local_id-1)*numfaces) == -DIRICHLET_BC).or. &
                 (e2n_local(iface + (local_id-1)*numfaces) == -HYDROSTATIC_BC))  then
            bound_id(iface) = 1
        end if
        if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
            neig_LP_id(iface) = cell_LP_id
        else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
                if (ghosted_id == conn%id_up(jface)) then
                    neig_LP_id(iface) = grid%ngmax_faces + conn%id_dn(jface) - 1
                else
                    neig_LP_id(iface) = grid%ngmax_faces + conn%id_up(jface) - 1
                end if
        end if
        ghosted_LP_id(numfaces+1+iface) = neig_LP_id(iface)

  end do
    ghosted_LP_id(numfaces + 1) = grid%ngmax_faces + ghosted_id - 1

!     write(*,*) ( ghosted_LP_id(iface),iface=1,7)  

    J = 0.
!    Jbl = 0.

   call MFDAuxJacobianLocal_LP( grid, aux_var, &
                                       rich_aux_vars(ghosted_id), global_aux_vars(ghosted_id), &
                                       sq_faces, option, J)
!    if (ghosted_id == 55) then
!    do iface=1,7 
!      write(*,*) (J(iface + (jface - 1)*7), jface = 1, 7)
!    end do
!    write(*,*)
!    end if

!    write(*,*) (bound_id(iface),iface=1,6)

   do jface = 1, numfaces 
        if (bound_id(jface) == 1) then
           do iface = 1, numfaces + 1
              J((iface-1)*nrow + jface) = 0.
              J((jface-1)*nrow + iface) = 0.
           end do
           J(jface + (jface - 1)*nrow) = 1.
        end if
   end do
 
! do jface = 1, numfaces + 1
!    do iface = 1, numfaces + 1
!       Jbl(iface + (jface-1)*(numfaces + 1)) = J(iface + (numfaces + 1)*(jface-1))
!    end do
!  end do
!
!   if (ghosted_id == 43)  then
!     do iface = 1, numfaces 
!      J(nrow + (iface-1)*nrow) = 0.
!     end do
!   end if


     call MatSetValuesLocal(A, numfaces + 1, ghosted_LP_id, numfaces + 1 , ghosted_LP_id, &
                                         J, ADD_VALUES,ierr)

!     call MatSetValuesLocal(A, numfaces + 1, ghosted_LP_id, numfaces + 1 + numfaces, ghosted_LP_id, &
!                                         Jbl, ADD_VALUES,ierr)


     call MatSetValuesLocal(A, 1, cell_LP_id, numfaces, neig_LP_id, aux_var%dRp_dneig, ADD_VALUES,ierr)

!     do iface = 1, numfaces
!     do iface = 1, 1
!        call MatSetValuesLocal(A, 1, cell_LP_id, 1, neig_LP_id(iface), aux_var%dRp_dneig(iface), ADD_VALUES,ierr)
!     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do



    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)


    deallocate(J)
!    deallocate(Jbl)
    deallocate(bound_id)
    deallocate(ghosted_LP_id)
    deallocate(neig_LP_id)
    deallocate(sq_faces)


  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_mfd.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

 

     call VecRestoreArrayF90(grid%e2n, e2n_local, ierr)

 


#endif

#if 0
  write(*,*) "EXIT MFDLP JACOBIAN"
  write(*,*) "richard 5932"
  read(*,*)
!   stop
#endif

end subroutine RichardsJacobianPatchMFDLP


! ************************************************************************** !
!
! RichardsCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsCreateZeroArray(patch,option)

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
!
! RichardsMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/15/08
!
! ************************************************************************** !
subroutine RichardsMaxChange(realization)

  use Realization_class
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization
  
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

#ifdef DASVYAT
!     call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx_faces.out',viewer,ierr)
!     call VecView(field%flow_dxx_faces, viewer,ierr)
!     call VecView(field%flow_dxx, viewer,ierr)
!     call PetscViewerDestroy(viewer, ierr)
!     write(*,*) "write flow_dxx_faces.out"
!     read(*,*)
#endif

  else 

     call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
     call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

  end if

end subroutine RichardsMaxChange

! ************************************************************************** !
!
! RichardsGetTecplotHeader: Returns Richards Lite contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
function RichardsGetTecplotHeader(realization,icolumn)
  
  use Realization_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: RichardsGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  option => realization%option
  field => realization%field
  
  string = ''
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')') 
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-sl"'')') icolumn
  else
    write(string2,'('',"sl"'')') 
  endif
  string = trim(string) // trim(string2)
 
  RichardsGetTecplotHeader = string

end function RichardsGetTecplotHeader

! ************************************************************************** !
!
! RichardsSetPlotVariables: Adds variables to be printed to list
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine RichardsSetPlotVariables(realization)
  
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
!
! RichardsPrintAuxVars: Prints out the contents of an auxvar
! author: Glenn Hammond
! date: 02/21/12
!
! ************************************************************************** !
subroutine RichardsPrintAuxVars(richards_auxvar,global_auxvar,cell_id)

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
!
! RichardsDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsDestroy(realization)

  use Realization_class

  type(realization_type) :: realization
  
  call RichardsDestroyPatch(realization)

end subroutine RichardsDestroy

! ************************************************************************** !
!
! RichardsDestroyPatch: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/03/09
!
! ************************************************************************** !
subroutine RichardsDestroyPatch(realization)

  use Realization_class

  implicit none

  type(realization_type) :: realization
  
  ! taken care of in auxiliary.F90

end subroutine RichardsDestroyPatch

end module Richards_module
