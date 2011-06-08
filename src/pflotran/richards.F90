module Richards_module

  use Richards_Aux_module
  use Global_Aux_module
#ifdef GLENN
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
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public RichardsResidual,RichardsJacobian, &
         RichardsUpdateFixedAccum,RichardsTimeCut,&
         RichardsSetup, RichardsNumericalJacTest, &
         RichardsInitializeTimestep, RichardsUpdateAuxVars, &
         RichardsMaxChange, RichardsUpdateSolution, &
         RichardsGetTecplotHeader, RichardsComputeMassBalance, &
         RichardsDestroy, RichardsResidualMFD, RichardsJacobianMFD, &
         RichardsInitialPressureReconstruction

contains

! ************************************************************************** !
!
! RichardsTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

  
  call VecCopy(field%flow_yy,field%flow_xx,ierr)

  if (option%mimetic) then
    call VecCopy(field%flow_yy_faces, field%flow_xx_faces, ierr)
    call RichardsUpdateAuxVars(realization)
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
      call RichardsSetupPatch(realization)
!      if (realization%option%mimetic) then
!         call RealizationSetUpBC4Faces(realization)
!      end if
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo


end subroutine RichardsSetup

! ************************************************************************** !
!
! RichardsSetupPatch: Creates arrays for auxilliary variables
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsSetupPatch(realization)

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
  PetscInt :: i, ierr
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)  
  
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
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  if (sum_connection > 0) then
    allocate(rich_aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RichardsAuxVarInit(rich_aux_vars_bc(iconn),option)
    enddo
    patch%aux%Richards%aux_vars_bc => rich_aux_vars_bc
  endif
  patch%aux%Richards%num_aux_bc = sum_connection
  
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call RichardsCreateZeroArray(patch,option)


end subroutine RichardsSetupPatch

! ************************************************************************** !
!
! RichardsComputeMassBalance: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsComputeMassBalance(realization,mass_balance)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
  
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
      call RichardsComputeMassBalancePatch(realization,mass_balance)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RichardsComputeMassBalance

! ************************************************************************** !
!
! RichardsComputeMassBalancePatch: Initializes mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RichardsComputeMassBalancePatch(realization,mass_balance)
 
  use Realization_module
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
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
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
  use Realization_module
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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
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

     write(*,*) "Mass conservation violation", max_violation
 
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
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

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

end subroutine RichardsZeroMassBalDeltaPatch

! ************************************************************************** !
!
! RichardsUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 12/19/08
!
! ************************************************************************** !
subroutine RichardsUpdateMassBalancePatch(realization)
 
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

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

end subroutine RichardsUpdateMassBalancePatch

! ************************************************************************** !
!
! RichardsUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateAuxVars(realization)

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
      if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
         if (.not.cur_patch%aux%Richards%aux_vars_cell_pressures_up_to_date) then
            call RichardsUpdateCellPressurePatch(realization)
         end if
      end if  
      call RichardsUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RichardsUpdateAuxVars


! ************************************************************************** !
!
! RichardsInitialPressureReconstruction: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



subroutine RichardsInitialPressureReconstruction(realization)

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
      call RichardsInitialPressureReconstructionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RichardsInitialPressureReconstruction



! ************************************************************************** !
!
! RichardsUpdateCellPressure: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



subroutine RichardsUpdateCellPressure(realization)

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
      if (.not.cur_patch%aux%Richards%aux_vars_cell_pressures_up_to_date) then
         call RichardsUpdateCellPressurePatch(realization)
      end if
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RichardsUpdateCellPressure 

! ************************************************************************** !
!
! RichardsUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsUpdateAuxVarsPatch(realization)

  use Realization_module
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
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)  
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)  
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase, i
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)  
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
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
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  



  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

   
    call RichardsAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id),rich_aux_vars(ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
  enddo

  call PetscLogEventEnd(logging%event_r_auxvars,ierr)

  call PetscLogEventBegin(logging%event_r_auxvars_bc,ierr)

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

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)  
  
  patch%aux%Richards%aux_vars_up_to_date = PETSC_TRUE

  call PetscLogEventEnd(logging%event_r_auxvars_bc,ierr)

end subroutine RichardsUpdateAuxVarsPatch

! ************************************************************************** !
!
! RichardsInitialPressureReconstructionPatch: Computes cell-centered pressures and 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !


subroutine RichardsInitialPressureReconstructionPatch(realization)


  use Realization_module
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
  type(discretization_type), pointer :: discretization

  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(mfd_auxvar_type), pointer :: mfd_aux_var

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
  PetscInt :: iphasebc, iphase
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_p(:), xx_loc_faces_p(:)
  PetscReal, pointer :: sq_faces(:), faces_pr(:)
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr



#ifdef DASVYAT_DEBUG
  write(*,*) "ENTER RichardsInitialPressureReconstructionPatch"
#endif


  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field


  rich_aux_vars => patch%aux%Richards%aux_vars
  global_aux_vars => patch%aux%Global%aux_vars

  call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces, ierr)

    
  call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)


  numfaces = 6     ! hex only
  allocate(sq_faces(numfaces))
  allocate(faces_pr(numfaces))


  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)   
 
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    mfd_aux_var => grid%MFD%aux_vars(local_id)
    do j = 1, mfd_aux_var%numfaces
       ghost_face_id = mfd_aux_var%face_id_gh(j)
       cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
       jface = grid%faces(ghost_face_id)%id
       sq_faces(j) = cur_connection_set%area(jface)

!       if (cur_connection_set%itype == INTERNAL_CONNECTION_TYPE) then
!         xx_loc_faces_p(ghost_face_id) = 1. !  DASVYAT test
!       end if 

       faces_pr(j) = xx_loc_faces_p(ghost_face_id)
    end do
 
        
        !geh - Ignore inactive cells with inactive materials
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif

        Res(1) = 0

        source_f = 0.
!        Res(1) = 0.

        call MFDAuxReconstruct(faces_pr, source_f, mfd_aux_var,&
                               rich_aux_vars(ghosted_id),global_aux_vars(ghosted_id), Res, &
                               sq_faces, option, xx_p(local_id:local_id))

   
!       write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)

  enddo


 patch%aux%Richards%aux_vars_cell_pressures_up_to_date = PETSC_TRUE 


  deallocate(sq_faces)
  deallocate(faces_pr)


  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p, ierr)

  call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF) 


!  call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_xx_loc_faces, field%flow_xx_faces, &   ! DASVYAT test
!                                INSERT_VALUES,SCATTER_REVERSE, ierr)
!   call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_xx_loc_faces, field%flow_xx_faces,&    ! DASVYAT test 
!                                INSERT_VALUES,SCATTER_REVERSE, ierr)


  
!   write(*,*) "DiscretizationGlobalToLocal"
!   read(*,*)

end subroutine RichardsInitialPressureReconstructionPatch

! ************************************************************************** !
!
! RichardsUpdateCellPressurePatch: Computes cell-centered pressures 
! author: Daniil Svyatskiy
! date: 12/01/10
!
! ************************************************************************** !



subroutine RichardsUpdateCellPressurePatch(realization)


  use Realization_module
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
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
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
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif
        

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
! RichardsUpdateAuxVarsPatchMFD: Computes cell-centered pressures and 
!                        updates the auxilliary variables associated with 
!                        the Richards problem
! author: Daniil Svyatskiy
! date: 07/29/10
!
! ************************************************************************** !

subroutine RichardsUpdateAuxVarsPatchMFD(realization)

  use Realization_module
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
  PetscInt :: iphasebc, iphase
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_loc_p(:), xx_faces_p(:), bc_loc_p(:), accum_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:), volume_p(:)  
  PetscReal, pointer :: sq_faces(:), darcy_v(:), faces_pr(:)
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr


! For Boundary Conditions
  PetscReal :: xxbc(realization%option%nflowdof)
  type(richards_auxvar_type), pointer :: rich_aux_vars_bc(:) 
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)


!#ifdef DASVYAT
  write(*,*) "ENTER RichardsUpdateAuxVarsPatchMFD"
!#endif

  call PetscLogEventBegin(logging%event_r_auxvars,ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field


  rich_aux_vars => patch%aux%Richards%aux_vars
  rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

    
  call GridVecGetArrayF90(grid,field%flow_xx, xx_loc_p, ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)
  call GridVecGetArrayF90(grid,field%perm0_xx, perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity0, porosity_loc_p,ierr)  
  call GridVecGetArrayF90(grid,field%volume, volume_p,ierr)  
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p,ierr)  


  numfaces = 6     ! hex only
  allocate(sq_faces(numfaces))
  allocate(faces_pr(numfaces))


  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)   
 
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    mfd_aux_var => grid%MFD%aux_vars(local_id)
    do j = 1, mfd_aux_var%numfaces
       ghost_face_id = mfd_aux_var%face_id_gh(j)
       cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
       jface = grid%faces(ghost_face_id)%id
       sq_faces(j) = cur_connection_set%area(jface)
       faces_pr(j) = xx_faces_p(ghost_face_id)
    end do
 
        
        !geh - Ignore inactive cells with inactive materials
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif
        call RichardsAccumulation(rich_aux_vars(ghosted_id), &
                                global_aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                volume_p(local_id), &
                                option,Res)

!        write(*,*) "accum_p", accum_p(local_id), "Res ",Res(1) , "diff", accum_p(local_id) - Res(1)

!       write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)
       
!       write(*,*) volume_p(local_id)*porosity_loc_p(ghosted_id)*global_aux_vars(ghosted_id)%den(1)/option%flow_dt 

        Res(1) = accum_p(local_id) - Res(1)

        source_f = 0.
!        Res(1) = 0.

         write(*,*) "press", xx_loc_p(local_id)
        call MFDAuxReconstruct(faces_pr, source_f, mfd_aux_var,&
                               rich_aux_vars(ghosted_id),global_aux_vars(ghosted_id), Res, &
                               sq_faces, option, xx_loc_p(local_id:local_id))


         write(*,*) "press", xx_loc_p(local_id)
   
       call RichardsAuxVarCompute(xx_loc_p(local_id:local_id),rich_aux_vars(ghosted_id), &
                       global_aux_vars(ghosted_id), &
                       patch%saturation_function_array((patch%sat_func_id(ghosted_id))%ptr, &
                       porosity_loc_p(local_id),perm_xx_loc_p(local_id), &                       
                       option)

!       write(*,*) "Sat", global_aux_vars(ghosted_id)%sat(1), "Pres", global_aux_vars(ghosted_id)%pres(1)

  enddo


  deallocate(sq_faces)
  deallocate(faces_pr)

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

      select case(boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
          xxbc(1) = boundary_condition%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          xxbc(1) = xx_loc_p(ghosted_id)
        case(UNIT_GRADIENT_BC)
          ! the auxvar is not needed for unit gradient
          cycle
      end select
      
      call RichardsAuxVarCompute(xxbc(1),rich_aux_vars_bc(sum_connection), &
                         global_aux_vars_bc(sum_connection), &
                         patch%saturation_function_array((patch%sat_func_id(ghosted_id))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                         
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call PetscLogEventEnd(logging%event_r_auxvars,ierr)

  call GridVecRestoreArrayF90(grid,field%flow_xx, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm0_xx, perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity0, porosity_loc_p,ierr)  
  call GridVecRestoreArrayF90(grid,field%volume, volume_p,ierr)  
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p,ierr)  

#ifdef DASVYAT_DEBUG
  write(*,*) "EXIT RichardsUpdateAuxVarsPatchMFD"
!  read(*,*)
#endif

#endif

end subroutine RichardsUpdateAuxVarsPatchMFD

! ************************************************************************** !
!
! RichardsInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine RichardsInitializeTimestep(realization)

  use Realization_module
  use MFD_module
  
  implicit none
  


  type(realization_type) :: realization

  PetscViewer :: viewer
  PetscErrorCode :: ierr

  call RichardsUpdateFixedAccum(realization)

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

  use Realization_module
  use Field_module
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization

  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
 
  if (realization%option%mimetic) then 
     call VecCopy(field%flow_xx_faces, field%flow_yy_faces, ierr)  
  end if

     call VecCopy(field%flow_xx,field%flow_yy,ierr)   

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RichardsUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  
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

  use Realization_module
    
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
      call RichardsUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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

  use Realization_module
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
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
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
  PetscReal :: porXvol

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert
  type(global_auxvar_type) :: global_aux_var_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  porXvol = por*vol/option%flow_dt
      
  ! accumulation term units = dkmol/dp
  J(1,1) = (global_aux_var%sat(1)*rich_aux_var%dden_dp+ &
            rich_aux_var%dsat_dp*global_aux_var%den(1))* &
           porXvol

  if (option%numerical_derivatives) then
    call GlobalAuxVarInit(global_aux_var_pert,option)  
    call RichardsAuxVarCopy(rich_aux_var,rich_aux_var_pert,option)
    call GlobalAuxVarCopy(global_aux_var,global_aux_var_pert,option)
    x(1) = global_aux_var%pres(1)
    call RichardsAccumulation(rich_aux_var,global_aux_var,por,vol,option,res)
    ideriv = 1
    pert = x(ideriv)*perturbation_tolerance
    x_pert = x
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
    call RichardsAccumulation(rich_aux_var_pert,global_aux_var,por,vol, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarDestroy(global_aux_var_pert)  
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
  PetscReal :: vol, por
       
  ! accumulation term units = kmol/s

#ifdef DASVYAT
!  write(*,*) "sat ",global_aux_var%sat(1), " den ", global_aux_var%den(1)
!  write(*,*) "por ", por, " volume ", vol, " dt ", option%flow_dt
#endif

    Res(1) = global_aux_var%sat(1) * global_aux_var%den(1) * por * vol / &
           option%flow_dt

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
                                  area, norm, dist_gravity,upweight, &
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
  PetscReal :: v_darcy, area, norm(3)
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
      if (dabs(norm(1))==1) then
        ukvr = rich_aux_var_up%kvr_x
        dukvr_dp_up = rich_aux_var_up%dkvr_x_dp
      else if (dabs(norm(2))==1) then
        ukvr = rich_aux_var_up%kvr_y
        dukvr_dp_up = rich_aux_var_up%dkvr_y_dp
      else if (dabs(norm(3))==1) then
        ukvr = rich_aux_var_up%kvr_z
        dukvr_dp_up = rich_aux_var_up%dkvr_z_dp
      end if
    else
      if (dabs(norm(1))==1) then
        ukvr = rich_aux_var_dn%kvr_x
        dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
      else if (dabs(norm(2))==1) then
        ukvr = rich_aux_var_dn%kvr_y
        dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
      else if (dabs(norm(3))==1) then
        ukvr = rich_aux_var_dn%kvr_z
        dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
      end if
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

  if (option%numerical_derivatives) then
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
                      area, norm, dist_gravity,upweight, &
                      option,v_darcy,res)
    ideriv = 1
    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_dn = x_dn(ideriv)*perturbation_tolerance
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
                      area, norm, dist_gravity,upweight, &
                      option,v_darcy,res_pert_up)
    call RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_pert_dn,global_aux_var_pert_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, norm, dist_gravity,upweight, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarDestroy(global_aux_var_pert_up)
    call GlobalAuxVarDestroy(global_aux_var_pert_dn)    
  endif

end subroutine RichardsFluxDerivative

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
                        area, norm, dist_gravity,upweight, &
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
  PetscReal :: v_darcy,area, norm(3)
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
  fluxm = 0.d0
  v_darcy = 0.D0  
  
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
!      ukvr = rich_aux_var_up%kvr
      if (dabs(dabs(norm(1))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_z
      end if
    else
!      ukvr = rich_aux_var_dn%kvr
      if (dabs(dabs(norm(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
    endif      
   

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi

!      write(*,*) "Gravity Input", Dq*ukvr*gravity
!	  write(*,*) "phi", global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)
!	  write(*,*) "Dq", Dq, "ukvr ", ukvr
   
      q = v_darcy * area

      fluxm = q*density_ave       
    endif
  endif 

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsFlux

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

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert_dn, rich_aux_var_pert_up
  type(global_auxvar_type) :: global_aux_var_pert_dn, global_aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  
  v_darcy = 0.d0
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
        
        density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
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
!          ukvr = rich_aux_var_up%kvr
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_x
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_y
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_z
          end if
        else
!          ukvr = rich_aux_var_dn%kvr
!          dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
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

        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
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
     
        if (ukvr*perm_dn>floweps) then
          v_darcy = perm_dn * ukvr * dphi
          q = v_darcy*area
          dq_dp_dn = perm_dn*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif 
     
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  if (option%numerical_derivatives) then
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
    pert_dn = x_dn(ideriv)*perturbation_tolerance    
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
    call GlobalAuxVarDestroy(global_aux_var_pert_up)
    call GlobalAuxVarDestroy(global_aux_var_pert_dn)      
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
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

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
!        ukvr = rich_aux_var_up%kvr
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_z
         end if
       else
!        ukvr = rich_aux_var_dn%kvr
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_z
         end if
       endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
!          write(*,*) "gravity ", gravity * Dq * ukvr
!          write(*,*) "phi", global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)
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
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(dist(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(dist(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
     
      if (ukvr*perm_dn>floweps) then
        v_darcy = perm_dn * ukvr * dphi
      endif 

  end select

#ifdef DASVYAT
!  write(*,*) "bound flux", v_darcy
#endif

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

  use Realization_module
  use Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module

  implicit none
  interface
     subroutine samrpetscobjectstateincrease(vec)
       implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       Vec :: vec
     end subroutine samrpetscobjectstateincrease
     
     subroutine SAMRCoarsenFaceFluxes(p_application, vec, ierr)
       implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: p_application
       Vec :: vec
       PetscErrorCode :: ierr
     end subroutine SAMRCoarsenFaceFluxes
  end interface

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
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
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RichardsResidualPatch1(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! now coarsen all face fluxes in case we are using SAMRAI to 
  ! ensure consistent fluxes at coarse-fine interfaces
  if(option%use_samr) then
     call SAMRCoarsenFaceFluxes(discretization%amrgrid%p_application, field%flow_face_fluxes, ierr)

     cur_level => realization%level_list%first
     do
        if (.not.associated(cur_level)) exit
        cur_patch => cur_level%patch_list%first
        do
           if (.not.associated(cur_patch)) exit
           realization%patch => cur_patch
           call RichardsResidualFluxContribPatch(r,realization,ierr)
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
      call RichardsResidualPatch2(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif
   
  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rresidual.bin',FILE_MODE_WRITE,&
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rxx.bin',FILE_MODE_WRITE,&
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  endif

call PetscLogEventEnd(logging%event_r_residual,ierr)


#ifdef DASVYAT_DEBUG
   call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
   call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  write(*,*) "End of RichardsResidual"
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

  use Realization_module
  use Field_module
  use Patch_module
  use Level_module
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
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: conn
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option

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
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

 ! ! now coarsen all face fluxes in case we are using SAMRAI to 
 ! ! ensure consistent fluxes at coarse-fine interfaces
 ! if(option%use_samr) then
 !    call SAMRCoarsenFaceFluxes(discretization%amrgrid%p_application, field%flow_face_fluxes, ierr)
!
!     cur_level => realization%level_list%first
!     do
!        if (.not.associated(cur_level)) exit
!        cur_patch => cur_level%patch_list%first
!        do
!           if (.not.associated(cur_patch)) exit
!           realization%patch => cur_patch
!           call RichardsResidualFluxContribPatch(r,realization,ierr)
!           cur_patch => cur_patch%next
!        enddo
!        cur_level => cur_level%next
!     enddo
!  endif


  ! pass #2 for boundary and source data
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RichardsResidualPatchMFD2(snes,xx,r,realization,ierr)
!      call RichardsCheckMassBalancePatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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

#ifdef DASVYAT_DEBUG
   write(*,*) "End RichardsResidualMFD"
   read(*,*) 
#endif   
!   write(*,*) "End RichardsResidualMFD"
!   read(*,*) 
!   stop

#endif
  
end subroutine RichardsResidualMFD

! ************************************************************************** !
!
! RichardsResidualFluxContribPatch: should be called only for SAMR
! author: Bobby Philip
! date: 02/17/09
!
! ************************************************************************** !
subroutine RichardsResidualFluxContribPatch(r,realization,ierr)
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
!!$
!!$  do axis=0,2  
!!$     call GridVecRestoreArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
!!$  enddo

end subroutine RichardsResidualFluxContribPatch

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
  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
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
  patch%aux%Richards%aux_vars_cell_pressures_up_to_date  = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
    call RichardsZeroMassBalDeltaPatch(realization)
  endif

#ifdef DASVYAT_DEBUG
      
!   write(*,*) "Internal faces"
!   do i = 1, grid%internal_connection_set_list%first%num_connections
!      write(*,*) "int_flux", option%myrank, i, patch%internal_velocities(option%nphase, i) 
!   end do  

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_xx_init.out', &
                              viewer,ierr)   
    call VecView(realization%field%flow_xx,viewer,ierr) 
!    call VecView(r,viewer,ierr) 

    call PetscViewerDestroy(viewer,ierr)

#endif



! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  !print *,' Finished scattering non deriv'

  if (option%use_samr) then
     do axis=0,2  
        call GridVecGetArrayF90(grid,axis,field%flow_face_fluxes, fluxes(axis)%flux_p, ierr)  
     enddo
  endif

  r_p = 0.d0
#if 1
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
                        cur_connection_set%area(iconn), cur_connection_set%dist(1:3,iconn), distance_gravity, &
                        upweight,option,v_darcy,Res)

      patch%internal_velocities(1,sum_connection) = v_darcy

#ifdef DASVYAT

!      write(*,*) "int flux ", sum_connection, v_darcy      
#endif

      
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
        fluxes(direction)%flux_p(flux_id) = Res(1)
      endif
      
#ifdef COMPUTE_INTERNAL_MASS_FLUX
      global_aux_vars(local_id_up)%mass_balance_delta(1,1) = &
        global_aux_vars(local_id_up)%mass_balance_delta(1,1) - Res(1)
#endif

      if(.not.option%use_samr) then
         
         if (local_id_up>0) then
            r_p(local_id_up) = r_p(local_id_up) + Res(1)
         endif
         
         if (local_id_dn>0) then
            r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
         endif
      endif
    enddo

    cur_connection_set => cur_connection_set%next
  enddo    

#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set


!    write(*,*) (boundary_condition%flow_aux_real_var(1, iconn), &
!              iconn =1,cur_connection_set%num_connections)
!    write(*,*)
     
    
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


#ifdef DASVYAT
!      write(*,*) "bound flux ", sum_connection, "fl",v_darcy, &
!            "lm", boundary_condition%flow_aux_real_var(1, iconn), &
!            "p", global_aux_vars(ghosted_id)%pres(1)
!            cur_connection_set%cntr(1, iconn), &      
!            cur_connection_set%cntr(2, iconn), &      
!            cur_connection_set%cntr(3, iconn)      
!            cur_connection_set%id_dn(iconn)       
#endif

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_aux_vars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_aux_vars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
        ! contribution to internal 
!        global_aux_vars(ghosted_id)%mass_balance_delta(1) = &
!          global_aux_vars(ghosted_id)%mass_balance_delta(1) + Res(1)
      endif

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
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(EAST_FACE)
              ghosted_id = ghosted_id+1
              flux_id = ((ghosted_id/ngxy)-1)*(nlx+1)*nly + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*(nlx+1)
              fluxes(direction)%flux_p(flux_id) = -Res(1)
           case(SOUTH_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(NORTH_FACE)
              ghosted_id = ghosted_id+ngx
              flux_id = ((ghosted_id/ngxy)-1)*nlx*(nly+1) + &
                        ((mod(ghosted_id,ngxy))/ngx-1)*nlx + &
                        mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -Res(1)
           case(BOTTOM_FACE)
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = Res(1)
           case(TOP_FACE)
              ghosted_id = ghosted_id+ngxy
              flux_id = ((ghosted_id/ngxy)-1)*nlx*nly &
                       +((mod(ghosted_id,ngxy))/ngx-1)*nlx &
                       +mod(mod(ghosted_id,ngxy),ngx)-1
              fluxes(direction)%flux_p(flux_id) = -Res(1)
         end select

!         fluxes(direction)%flux_p(flux_id) = Res(1)

      else
         r_p(local_id)= r_p(local_id) - Res(1)
      endif

   enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  

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
  type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  
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
  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  ! Accumulation terms ------------------------------------
  if (.not.option%steady_state) then
#if 1
    

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
#endif
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
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
      end select
!      if (option%compute_mass_balance_new) then
        ! need to added global aux_var for src/sink
!        global_aux_vars_ss(ghosted_id)%mass_balance_delta(1) = &
!          global_aux_vars_ss(ghosted_id)%mass_balance_delta(1) - qsrc_kg
!      endif
      r_p(local_id) = r_p(local_id) - qsrc_mol
    enddo
    source_sink => source_sink%next
  enddo
#endif


!#ifdef DASVYAT
!     write(*,*) "richards 2608"
!     stop
!#endif

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

  use Connection_module
  use Realization_module
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
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3), den, ukvr
  PetscInt :: icell, iface, jface, i,j, numfaces, ghost_face_id, ghost_face_jd
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
!  allocate(Smatrix(numfaces, numfaces))

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
    end do

    ghosted_id = grid%nL2G(icell)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif


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

      call MFDAuxFluxes(patch, grid, ghosted_id, xx_p(icell:icell), face_pr, aux_var, PermTensor, &
                        rich_aux_vars(ghosted_id), global_aux_vars(ghosted_id), &
                        sq_faces, option)

      
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

#ifdef DASVYAT
!  write(*,*) "richards 2822"
!  write(*,*) "End RichardsResidualPatchMFD1"
!  stop
#endif

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

  use Connection_module
  use Realization_module
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


  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:),  accum_p(:), bc_faces_p(:)
  PetscReal, pointer ::  perm_xx_loc_p(:),  perm_yy_loc_p(:), perm_zz_loc_p(:), flow_xx_p(:)
  PetscReal, pointer ::  perm_xz_loc_p(:),  perm_xy_loc_p(:), perm_yz_loc_p(:)
  PetscReal, pointer :: xx_loc_faces_p(:)


  PetscReal :: qsrc, qsrc_mol


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
  PetscReal :: sq_faces(6), rhs(6), bc_g(6), bc_h(6), face_pres(6), bnd(6)
  PetscReal :: Accum(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3)
  PetscInt :: icell, iface, jface, numfaces, ghost_face_id, ghost_face_jd
  PetscScalar, pointer :: e2n_local(:)
  


  
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
  if(ierr /= 0) write(*,*)"sq_faces allocation error" 

  allocate(face_pres(numfaces), stat = ierr)
  if(ierr /= 0) write(*,*)"sq_faces allocation error" 

  allocate(rhs(numfaces), stat = ierr)
  if(ierr /= 0) write(*,*)"rhs allocation error" 

  allocate(bc_g(numfaces), stat = ierr)
  if(ierr /= 0) write(*,*)"bc_g allocation error" 

  allocate(bc_h(numfaces), stat = ierr)
  if(ierr /= 0) write(*,*)"bc_h allocation error" 

  allocate(bnd (numfaces), stat = ierr)
  if(ierr /= 0) write(*,*)"bnd allocation error" 
#endif
  stride = 6 !hex only


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

        end do
         

        !geh - Ignore inactive cells with inactive materials
        if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
        endif
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

 

        call MFDAuxGenerateRhs(grid, ghosted_id, PermTensor, bc_g, source_f, bc_h, aux_var, &
                                 rich_aux_vars(ghosted_id),&
                                 global_aux_vars(ghosted_id),&
                                 Accum, &
                                 porosity_loc_p(ghosted_id), volume_p(local_id),&
                                 flow_xx_p(local_id:local_id), face_pres, bnd,&                                 
                                 sq_faces, option, rhs) 

 
     
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
!
!		write(*,*) (rhs(iface),iface=1,6)
!		write(*,*)
        
  enddo

!  do iface =1,grid%ngmax_faces
!    write(*,*) "residual_p ", grid%fG2P(iface), r_p(iface), xx_loc_faces_p(iface)
!  end do

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

#endif


end subroutine RichardsResidualPatchMFD2


! ************************************************************************** !
!
! RichardsJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module

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
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif

      call RichardsJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
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
      if(option%use_samr) then
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif

      call RichardsJacobianPatch2(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

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
  
end subroutine RichardsJacobian
! ************************************************************************** !
!
! RichardsJacobian: Computes the Jacobian for MFD Discretization
! author: Daniil Svyatskiy
! date: 09/17/10
!
! ************************************************************************** !
subroutine RichardsJacobianMFD(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module

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
!      if(option%use_samr) then
!         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
!      endif

      call RichardsJacobianPatchMFD(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo


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
  use Realization_module
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

#ifdef GLENN
  if (option%use_matrix_buffer) then
    if (associated(patch%aux%Richards%matrix_buffer)) then
      call MatrixBufferZero(patch%aux%Richards%matrix_buffer)
    else
      patch%aux%Richards%matrix_buffer => MatrixBufferCreate()
      call MatrixBufferInit(A,patch%aux%Richards%matrix_buffer,grid)
    endif
  endif
#endif

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
                                  cur_connection_set%area(iconn), cur_connection_set%dist(1:3,iconn), distance_gravity, &
                                  upweight,option,&
                                  patch%saturation_function_array(icap_up)%ptr,&
                                  patch%saturation_function_array(icap_dn)%ptr,&
                                  Jup,Jdn)
      if (local_id_up > 0) then
#ifdef GLENN
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
#ifdef GLENN
        endif
#endif
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
#ifdef GLENN
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
#ifdef GLENN
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
#ifdef GLENN
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jdn(1,1))
      else
#endif
        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                               ADD_VALUES,ierr)
#ifdef GLENN
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
  use Realization_module
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  interface
     subroutine SAMRSetJacobianSourceOnPatch(which_pc, index, val, p_application, p_patch) 
#include "finclude/petscsysdef.h"

       PetscInt :: which_pc
       PetscInt :: index
       PetscReal :: val
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSourceOnPatch

     subroutine SAMRSetJacobianSrcCoeffsOnPatch(which_pc, p_application, p_patch) 
#include "finclude/petscsysdef.h"

       PetscInt :: which_pc
       PetscFortranAddr :: p_application
       PetscFortranAddr :: p_patch
     end subroutine SAMRSetJacobianSrcCoeffsOnPatch
  end interface

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
#ifdef GLENN
    if (option%use_matrix_buffer) then
      call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                           ghosted_id,Jup(1,1))
    else
#endif
      call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                             ADD_VALUES,ierr)
#ifdef GLENN
    endif
#endif
!!$    if(option%use_samr) then
!!$       flow_pc = 0
!!$       call SAMRSetJacobianSourceOnPatch(flow_pc, ghosted_id-1, Jup(1,1), &
!!$       realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
!!$    endif
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
    
    qsrc = source_sink%flow_condition%rate%dataset%cur_value(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
      
      Jup = 0.d0
      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS,SCALED_MASS_RATE_SS)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_aux_vars(ghosted_id)%dden_dp*FMWH2O
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          Jup(1,1) = -qsrc*rich_aux_vars(ghosted_id)%dden_dp*FMWH2O* &
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)

      end select
#ifdef GLENN
      if (option%use_matrix_buffer) then
        call MatrixBufferAdd(patch%aux%Richards%matrix_buffer,ghosted_id, &
                             ghosted_id,Jup(1,1))
      else
#endif
        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)  
#ifdef GLENN
      endif
#endif

!!$      if(option%use_samr) then
!!$         flow_pc = 0
!!$         call SAMRSetJacobianSourceOnPatch(flow_pc, ghosted_id-1, Jup(1,1), &
!!$         realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
!!$      endif
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

#ifdef GLENN
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
#ifdef GLENN
  if (.not.option%use_matrix_buffer) then
#endif
    if (patch%aux%Richards%inactive_cells_exist) then
      qsrc = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%Richards%n_zero_rows, &
                            patch%aux%Richards%zero_rows_local_ghosted, &
                            qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
    endif
#ifdef GLENN
  endif
#endif

  if(option%use_samr) then
     flow_pc = 0
     call SAMRSetJacobianSrcCoeffsOnPatch(flow_pc, &
          realization%discretization%amrgrid%p_application, grid%structured_grid%p_samr_patch)
  endif

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
  use Realization_module
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
     ukvr = rich_aux_vars(ghosted_id)%kvr_x 
     dkvr_x_dp = rich_aux_vars(ghosted_id)%dkvr_x_dp
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

#ifdef DASVYAT_DEBUG
  write(*,*) "EXIT MFD JACOBIAN"
  write(*,*) "richard 4587"
!  stop
#endif

end subroutine RichardsJacobianPatchMFD



! ************************************************************************** !
!
! RichardsCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsCreateZeroArray(patch,option)

  use Realization_module
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

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%nflowdof
      endif
    enddo
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
      endif
    enddo
  endif

  patch%aux%Richards%zero_rows_local => zero_rows_local
  patch%aux%Richards%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%Richards%n_zero_rows = n_zero_rows
  
  if(.not. (option%use_samr)) then
     call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                        MPI_MAX,option%mycomm,ierr)
     if (flag > 0) patch%aux%Richards%inactive_cells_exist = PETSC_TRUE
     
     if (ncount /= n_zero_rows) then
        print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
        stop
     endif
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

  use Realization_module
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
!     call VecStrideNorm(field%flow_dxx_faces,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

     call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
     call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

#ifdef DASVYAT_DEBUG
     call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx_faces.out',viewer,ierr)
     call VecView(field%flow_dxx_faces, viewer,ierr)
     call VecView(field%flow_dxx, viewer,ierr)
     call PetscViewerDestroy(viewer, ierr)
     write(*,*) "write flow_dxx_faces.out"
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
  
  use Realization_module
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
! RichardsDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsDestroy(realization)

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
      call RichardsDestroyPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RichardsDestroy

! ************************************************************************** !
!
! RichardsDestroyPatch: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/03/09
!
! ************************************************************************** !
subroutine RichardsDestroyPatch(realization)

  use Realization_module

  implicit none

  type(realization_type) :: realization
  
  ! taken care of in auxilliary.F90

end subroutine RichardsDestroyPatch

end module Richards_module
