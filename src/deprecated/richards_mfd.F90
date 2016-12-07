module Richards_MFD_module
#include "finclude/petscmat.h"
  use petscmat
  use Richards_Aux_module
  use Richards_Common_module
  use Global_Aux_module
#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsResidualMFD, RichardsJacobianMFD, &
         RichardsResidualMFDLP, RichardsJacobianMFDLP, &
         RichardsUpdateAuxVarsPatchMFDLP

contains

! ************************************************************************** !

subroutine RichardsCheckMassBalancePatch(realization)
#include "finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use MFD_Aux_module
  use MFD_module
  use Auxiliary_module
  use Material_Aux_class
  
  implicit none
  PetscInt :: local_id, ghosted_cell_id

  PetscReal, pointer :: accum_p(:), xx_faces_loc_p(:)

  type(realization_type) :: realization


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal :: mass_conserv, res(1)
  PetscErrorCode :: ierr

  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn
  PetscScalar :: sq_faces(6), den(6), dden_dp(6), max_violation
  PetscInt :: jface, j, numfaces, ghost_face_id, local_face_id, dir, bound_id
  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  field => realization%field

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_loc_p,  &
                      ierr);CHKERRQ(ierr)

  max_violation = 0.

  do local_id = 1, grid%nlmax
    if ((local_id.ne.1).and.(local_id.ne.150)) cycle

    mass_conserv = 0.
    auxvar => grid%MFD%auxvars(local_id)
    ghosted_cell_id = grid%nL2G(local_id)

    do j = 1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(j)
      local_face_id = grid%fG2L(ghost_face_id)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      jface = grid%faces(ghost_face_id)%id
      sq_faces(j) = conn%area(jface)
      call MFDComputeDensity(global_auxvars(ghosted_cell_id), xx_faces_loc_p(ghost_face_id), &
                                            den(j), dden_dp(j), option)
      dir = 1
      if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(jface)==ghosted_cell_id) dir = -1

      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        if (local_face_id > 0) then
          bound_id = grid%fL2B(local_face_id)
          if (bound_id>0) then
            mass_conserv = mass_conserv - den(j)*patch%boundary_velocities(option%nphase, bound_id)*sq_faces(j)
          endif
        endif
      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        mass_conserv = mass_conserv + den(j)*patch%internal_velocities(option%nphase, jface)*sq_faces(j)* dir
      endif
    enddo
        
    call RichardsAccumulation(rich_auxvars(ghosted_cell_id), &
                                global_auxvars(ghosted_cell_id), &
                                material_auxvars(ghosted_cell_id), &
                                option, res)

    mass_conserv = mass_conserv + res(1) - accum_p(local_id)
    max_violation = max(max_violation, abs(mass_conserv))

  enddo

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, xx_faces_loc_p,  &
                          ierr);CHKERRQ(ierr)

end subroutine RichardsCheckMassBalancePatch

! ************************************************************************** !

subroutine RichardsUpdateCellPressure(realization)
  ! 
  ! Computes cell-centered pressures
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 12/01/10
  ! 

  use Realization_class

  type(realization_type) :: realization

  if (.not.realization%patch%aux%Richards% &
      auxvars_cell_pressures_up_to_date) then
    call RichardsUpdateCellPressurePatch(realization)
  endif

end subroutine RichardsUpdateCellPressure 

! ! ************************************************************************** !
! 
! subroutine RichardsInitialPressureReconstructionPatch(realization)
! !
! ! Computes cell-centered pressures
! ! Author: Daniil Svyatskiy
! ! Date: 12/01/10
! !
!   use Realization_class
!   use Discretization_module
!   use Patch_module
!   use Option_module
!   use Field_module
!   use Grid_module
!   use Connection_module
!   use MFD_module
!   use MFD_Aux_module
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
!   type(richards_auxvar_type), pointer :: rich_auxvars(:)
!   type(global_auxvar_type), pointer :: global_auxvars(:)
!   type(mfd_auxvar_type), pointer :: mfd_auxvar
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
!   rich_auxvars => patch%aux%Richards%auxvars
!   global_auxvars => patch%aux%Global%auxvars
! 
!   call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces, ierr)
! 
!     
!   call VecGetArrayF90(field%flow_xx, xx_p, ierr)
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
!     mfd_auxvar => grid%MFD%auxvars(local_id)
!     do j = 1, mfd_auxvar%numfaces
!        ghost_face_id = mfd_auxvar%face_id_gh(j)
!        cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
!        jface = grid%faces(ghost_face_id)%id
!        sq_faces(j) = cur_connection_set%area(jface)
! 
! !       if (cur_connection_set%itype == INTERNAL_CONNECTION_TYPE) then
! !         xx_loc_faces_p(ghost_face_id) = 1. !  DASVYAT test
! !       endif 
! 
!        faces_pr(j) = xx_loc_faces_p(ghost_face_id)
!     enddo
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
!         call MFDAuxReconstruct(faces_pr, source_f, mfd_auxvar,&
!                                rich_auxvars(ghosted_id),global_auxvars(ghosted_id), Res, &
!                                sq_faces, option, xx_p(local_id:local_id))
! 
!    
! !       write(*,*) "Sat", global_auxvars(ghosted_id)%sat(1), "Pres", global_auxvars(ghosted_id)%pres(1)
! 
!   enddo
! 
! 
!  patch%aux%Richards%auxvars_cell_pressures_up_to_date = PETSC_TRUE 
! 
! 
!   deallocate(sq_faces)
!   deallocate(faces_pr)
! 
! 
!   call VecRestoreArrayF90(field%flow_xx, xx_p, ierr)
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

subroutine RichardsUpdateCellPressurePatch(realization)
  ! 
  ! Computes cell-centered pressures
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 12/01/10
  ! 

  use Realization_class
  use Discretization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use MFD_module
  use MFD_Aux_module
  use Material_Aux_class
  
  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(mfd_auxvar_type), pointer :: mfd_auxvar

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
  PetscInt :: iphasebc, iphase
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_p(:), xx_loc_faces_p(:), work_loc_faces_p(:)
  PetscReal, pointer ::  accum_p(:)
  PetscReal, pointer :: sq_faces(:), faces_DELTA_pr(:), faces_pr(:) 
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  discretization => realization%discretization
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%work_loc_faces, work_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p,ierr);CHKERRQ(ierr)

  numfaces = 6     ! hex only
  allocate(sq_faces(numfaces))
  allocate(faces_DELTA_pr(numfaces))
  allocate(faces_pr(numfaces))

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)
 
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    mfd_auxvar => grid%MFD%auxvars(local_id)
    do j = 1, mfd_auxvar%numfaces
      ghost_face_id = mfd_auxvar%face_id_gh(j)
      cur_connection_set => grid%faces(ghost_face_id)%conn_set_ptr
      jface = grid%faces(ghost_face_id)%id
      sq_faces(j) = cur_connection_set%area(jface)
      faces_pr(j) = work_loc_faces_p(ghost_face_id)
       
      faces_pr(j) = xx_loc_faces_p(ghost_face_id)
       
      faces_DELTA_pr(j) = xx_loc_faces_p(ghost_face_id) - work_loc_faces_p(ghost_face_id)
    enddo

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
        

    call RichardsAccumulation(rich_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              option,Res)

    Res(1) = Res(1) - accum_p(local_id)

    source_f = 0.

    call MFDAuxUpdateCellPressure(faces_pr, faces_DELTA_pr, mfd_auxvar,&
                               option, xx_p(local_id:local_id))
  enddo
  patch%aux%Richards%auxvars_cell_pressures_up_to_date = PETSC_TRUE

  deallocate(sq_faces)
  deallocate(faces_DELTA_pr)
  deallocate(faces_pr)

  call DiscretizationGlobalToLocal(discretization, field%flow_xx, &
                                   field%flow_xx_loc, NFLOWDOF)
  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%work_loc_faces, work_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p,ierr);CHKERRQ(ierr)

end subroutine RichardsUpdateCellPressurePatch

! ************************************************************************** !

subroutine RichardsUpdateAuxVarsPatchMFDLP(realization)
  ! 
  ! Computes  updates the auxiliary variables associated with
  ! the Richards problem for LP formulation
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 07/29/10
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
  use MFD_module
  use MFD_Aux_module
  use Material_Aux_class
  
  implicit none

  type(realization_type) :: realization

#ifdef DASVYAT
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(mfd_auxvar_type), pointer :: mfd_auxvar
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, i,j
  PetscInt :: iphasebc, iphase, LP_cell_id, LP_face_id
  PetscInt :: ghost_face_id, iface, jface, numfaces
  PetscReal, pointer :: xx_p(:), xx_LP_loc_p(:), bc_loc_p(:)
  PetscReal, pointer :: sq_faces(:), darcy_v(:), faces_pr(:)
  PetscReal :: Res(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscErrorCode :: ierr

! For Boundary Conditions
  PetscReal :: xxbc(realization%option%nflowdof)
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  call PetscLogEventBegin(logging%event_r_auxvars,ierr);CHKERRQ(ierr)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  
  call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_LP_loc_p,  &
                      ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

     !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

 
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    LP_cell_id = grid%ngmax_faces + ghosted_id
    call RichardsAuxVarCompute(xx_LP_loc_p(LP_cell_id:LP_cell_id),rich_auxvars(ghosted_id), &
                       global_auxvars(ghosted_id), &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       material_auxvars(ghosted_id)%porosity, &
                       material_auxvars(ghosted_id)%permeability(perm_xx_index), &                       
                       option)

    local_id = grid%nG2L(ghosted_id)
    if (local_id > 0) xx_p(local_id) = xx_LP_loc_p(LP_cell_id)
  enddo

  boundary_condition => patch%boundary_condition_list%first
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
      
      call RichardsAuxVarCompute(xxbc(1),rich_auxvars_bc(sum_connection), &
                         global_auxvars_bc(sum_connection), &
                         patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                         material_auxvars(ghosted_id)%porosity, &
                         material_auxvars(ghosted_id)%permeability(perm_xx_index), &                       
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call PetscLogEventEnd(logging%event_r_auxvars,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_LP_loc_p,  &
                          ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsUpdateAuxVarsPatchMFDLP

! ************************************************************************** !

subroutine RichardsResidualMFD(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation for MFD discretization
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 05/26/10
  ! 

#include "finclude/petscsnes.h"
  use petscssnes
  use Logging_module
  use Realization_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Connection_module
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, PERMEABILITY_Z, &
                               PERMEABILITY_XY, PERMEABILITY_XZ, PERMEABILITY_YZ

  implicit none

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


  call PetscLogEventBegin(logging%event_r_residual,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  r_p = 0.
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)

  call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r, &
                                INSERT_VALUES,SCATTER_REVERSE,  &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r,&
                                INSERT_VALUES,SCATTER_REVERSE, ierr)


  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
   call DiscretizationGlobalToLocalFaces(discretization, xx, field%flow_xx_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocalFaces(discretization, r, field%flow_r_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF)

  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)

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
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_XZ,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_XZ,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_XY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_XY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_YZ,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_YZ,ZERO_INTEGER)

  if (realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
         call RichardsUpdateCellPressure(realization)
  endif  

  ! pass #1 for flux terms and accumulation term
  call RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)

  ! pass #2 for boundary and source data
  call RichardsResidualPatchMFD2(snes,xx,r,realization,ierr)

  call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces,  &
               ierr);CHKERRQ(ierr)

  call VecScatterBegin( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r, &
                                ADD_VALUES,SCATTER_REVERSE,  &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd ( discretization%MFD%scatter_gtol_faces, field%flow_r_loc_faces, r,&
                                ADD_VALUES,SCATTER_REVERSE, ierr)


  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RresidualMFD.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RxxMFD.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
 endif

  call PetscLogEventEnd(logging%event_r_residual,ierr);CHKERRQ(ierr)

#endif
  
end subroutine RichardsResidualMFD

! ************************************************************************** !

subroutine RichardsResidualMFDLP(snes,xx,r,realization,ierr)

#include "finclude/petscsnes.h"
  use petscsnes
  use Logging_module
  use Realization_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Connection_module
  use Mass_Transfer_module, only : mass_transfer_type  

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscReal, pointer :: v_p(:)
  PetscViewer :: viewer
  PetscErrorCode :: ierr

#if DASVYAT

  PetscInt :: i, iface, icell
  PetscReal, pointer :: flow_xx_loc_p(:), r_p(:)
  PetscReal :: rnorm
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: conn
  type(mass_transfer_type), pointer :: cur_mass_transfer
  
  field => realization%field
  discretization => realization%discretization
  option => realization%option


  call PetscLogEventBegin(logging%event_r_residual,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  r_p = 0.
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)

  call VecScatterBegin( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r, &
                                INSERT_VALUES,SCATTER_REVERSE,  &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd ( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r,&
                                INSERT_VALUES,SCATTER_REVERSE, ierr)

 
  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
   call DiscretizationGlobalToLocalLP(discretization, xx, field%flow_xx_loc_faces, NFLOWDOF)
   call DiscretizationGlobalToLocalLP(discretization, r, field%flow_r_loc_faces, NFLOWDOF)
!    call DiscretizationGlobalToLocal(discretization, field%flow_xx, field%flow_xx_loc, NFLOWDOF)

  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xz_loc,field%perm_xz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_xy_loc,field%perm_xy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yz_loc,field%perm_yz_loc,ONEDOF)

  ! pass #1 for flux terms and accumulation term
   call RichardsResidualPatchMFDLP1(snes,xx,r,realization,ierr)

  ! pass #2 for boundary and source data
   call RichardsResidualPatchMFDLP2(snes,xx,r,realization,ierr)
  !call RichardsCheckMassBalancePatch(realization)

  ! Mass Transfer
  if (associated(realization%flow_mass_transfer_list)) then
    cur_mass_transfer => realization%flow_mass_transfer_list
    call VecGetArrayF90(field%flow_r_loc_faces,r_p,ierr);CHKERRQ(ierr)
    do
      if (.not.associated(cur_mass_transfer)) exit
      call VecGetArrayF90(cur_mass_transfer%vec,v_p,ierr);CHKERRQ(ierr)
      do icell=1,realization%patch%grid%nlmax
        r_p(icell+realization%patch%grid%ngmax_faces) = &
            r_p(icell+realization%patch%grid%ngmax_faces) + v_p(icell)
      enddo
      call VecRestoreArrayF90(cur_mass_transfer%vec,v_p,ierr);CHKERRQ(ierr)
      cur_mass_transfer => cur_mass_transfer%next
    enddo
    call VecRestoreArrayF90(field%flow_r_loc_faces,r_p,ierr);CHKERRQ(ierr)
  endif  

   call VecScatterBegin( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r, &
                                ADD_VALUES,SCATTER_REVERSE,  &
                        ierr);CHKERRQ(ierr)
   call VecCopy(field%flow_xx_loc_faces, field%work_loc_faces,  &
                ierr);CHKERRQ(ierr)
   call VecScatterEnd ( discretization%MFD%scatter_gtol_LP, field%flow_r_loc_faces, r,&
                                ADD_VALUES,SCATTER_REVERSE, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RresidualMFD.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RxxMFD.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
 endif

  call PetscLogEventEnd(logging%event_r_residual,ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsResidualMFDLP

! ************************************************************************** !

subroutine RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)
  ! 
  ! Computes flux and accumulation
  ! terms of the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 
#include "finclude/petscvec.h"
  use petscvec
  
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
  use Material_Aux_class
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), xx_loc_faces_p(:), xx_p(:), &
                        work_loc_faces_p(:)

  PetscReal, pointer :: face_fluxes_p(:), bc_faces_p(:)

  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn, bound_id 

  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn
  PetscScalar, pointer :: sq_faces(:), e2n_local(:), Smatrix(:,:), face_pr(:)
  PetscScalar, pointer :: neig_den(:), neig_kvr(:), neig_pres(:), bnd(:)
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3), den, ukvr
  PetscInt :: icell, iface, jface, i,j, numfaces, ghost_face_id, ghost_face_jd
  PetscInt :: ghost_neig_id
  PetscViewer :: viewer

  call PetscLogEventBegin(logging%event_flow_residual_mfd1,  &
                          ierr);CHKERRQ(ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  discretization => realization%discretization 
  material_auxvars => patch%aux%Material%auxvars
  
  call RichardsUpdateAuxVarsPatchMFDLP(realization)

  patch%aux%Richards%auxvars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  patch%aux%Richards%auxvars_cell_pressures_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
!    call RichardsZeroMassBalDeltaPatch(realization)
  endif


! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr);CHKERRQ(ierr)

  numfaces = 6 ! hex only
  allocate(sq_faces(numfaces))
  allocate(face_pr(numfaces))
  allocate(neig_den(numfaces))
  allocate(neig_kvr(numfaces))
  allocate(neig_pres(numfaces))
  allocate(bnd(numfaces))

  bnd = 0

  do icell = 1, grid%nlmax

    auxvar => grid%MFD%auxvars(icell)
    do j = 1, auxvar%numfaces
      ghosted_id = grid%nL2G(icell)
      ghost_face_id = auxvar%face_id_gh(j)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      jface = grid%faces(ghost_face_id)%id
      sq_faces(j) = conn%area(jface)
      face_pr(j) = xx_loc_faces_p(ghost_face_id)
 
      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        neig_den(j) = global_auxvars(ghosted_id)%den(1)
      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        if (ghosted_id == conn%id_up(jface)) then
          ghost_neig_id = conn%id_dn(jface)
        else
          ghost_neig_id = conn%id_up(jface)
        endif
        neig_den(j) = global_auxvars(ghost_neig_id)%den(1)
      endif
    enddo

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_id)%permeability(perm_yz_index)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(3,2) = PermTensor(2,3)


    call PetscLogEventBegin(logging%event_flow_flux_mfd, ierr);CHKERRQ(ierr)
    call MFDAuxFluxes(patch, grid, ghosted_id, xx_p(icell:icell), face_pr, auxvar, PermTensor, &
                      rich_auxvars(ghosted_id), global_auxvars(ghosted_id), &
                      sq_faces,  bnd, neig_den, neig_kvr, neig_pres, option)
    call PetscLogEventEnd(logging%event_flow_flux_mfd, ierr);CHKERRQ(ierr)
      
  enddo

  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)

  deallocate(sq_faces)
  deallocate(face_pr)
  deallocate(neig_den)
  deallocate(bnd)

  call PetscLogEventEnd(logging%event_flow_residual_mfd1, ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsResidualPatchMFD1

! ************************************************************************** !

subroutine RichardsResidualPatchMFD2(snes,xx,r,realization,ierr)
  ! 
  ! Computes the boundary and source/sink terms of
  ! the residual equation for MFD discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

#include "finclude/petscvec.h"
  use petscvec
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
  use MFD_Aux_module
  use Material_Aux_class
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: i,j
  PetscInt :: local_id, ghosted_id
  
  PetscScalar, pointer :: r_p(:), accum_p(:), bc_faces_p(:)
  PetscScalar, pointer :: xx_loc_faces_p(:), flow_xx_p(:)

  PetscScalar :: qsrc, qsrc_mol

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink, boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection, bc_type, stride


  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn
! PetscReal, pointer :: sq_faces(:), rhs(:), bc_g(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: sq_faces(6), rhs(6), bc_g(6), bc_h(6), face_pres(6), bnd(6)
  PetscScalar :: neig_den(6), neig_kvr(6), neig_dkvr_dp(6)
  PetscScalar :: Accum(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscScalar :: Res(realization%option%nflowdof), PermTensor(3,3)
  PetscInt :: icell, iface, jface, numfaces, ghost_face_id, ghost_face_jd, ghost_neig_id
  PetscScalar, pointer :: e2n_local(:)
  


  call PetscLogEventBegin(logging%event_flow_residual_mfd2,  &
                          ierr);CHKERRQ(ierr)
  
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
  
! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

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

    auxvar => grid%MFD%auxvars(local_id)
    do j = 1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(j)
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
        bc_h(j) = bc_faces_p(ghost_face_id)
      else if (e2n_local(j + (local_id-1)*stride) < 0) then
        call printMsg(option,'Type of BC is not supported by MFD mode')
        stop
      endif

      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        neig_den(j) = global_auxvars_bc(jface)%den(1)
#ifdef USE_ANISOTROPIC_MOBILITY
        neig_kvr(j) = rich_auxvars(ghosted_id)%kvr_x
        neig_dkvr_dp(j) = rich_auxvars(ghosted_id)%dkvr_x_dp
#else
        neig_kvr(j) = rich_auxvars_bc(jface)%kvr
        neig_dkvr_dp(j) = 0.
#endif
      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        if (ghosted_id == conn%id_up(jface)) then
          ghost_neig_id = conn%id_dn(jface)
        else
          ghost_neig_id = conn%id_up(jface)
        endif
 
        neig_den(j) = global_auxvars(ghost_neig_id)%den(1)
#ifdef USE_ANISOTROPIC_MOBILITY
        neig_kvr(j) = rich_auxvars(ghost_neig_id)%kvr_x
        neig_dkvr_dp(j) = rich_auxvars(ghost_neig_id)%dkvr_x_dp
#else
        neig_kvr(j) = rich_auxvars(ghost_neig_id)%kvr
        neig_dkvr_dp(j) = rich_auxvars(ghost_neig_id)%dkvr_dp
#endif
      endif
    enddo

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsAccumulation(rich_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auvars(ghosted_id)%porosity, &
                              material_auvars(ghosted_id)%volume, &
                              option, Res)

    Accum(1) = Res(1) - accum_p(local_id)
 
    source_f = 0.
    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_id)%permeability(perm_yz_index)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,2) = PermTensor(2,3)

    call PetscLogEventBegin(logging%event_flow_rhs_mfd, ierr);CHKERRQ(ierr)
!    call MFDAuxGenerateRhs(patch, grid, ghosted_id, PermTensor, bc_g, source_f, bc_h, auxvar, &
!                                 rich_auxvars(ghosted_id),&
!                                 global_auxvars(ghosted_id),&
!                                 Accum, &
!                                 porosity_loc_p(ghosted_id), volume_p(local_id),&
!                                 flow_xx_p(local_id:local_id), face_pres, bnd,&                                 
!                                 sq_faces, neig_den, neig_kvr, neig_dkvr_dp, option, rhs) 
    call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr);CHKERRQ(ierr)
    call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr);CHKERRQ(ierr)

    do iface=1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(iface)

      if (e2n_local((local_id -1)*numfaces + iface) < 0) then
        if ((e2n_local((local_id -1)*numfaces + iface) == -DIRICHLET_BC).or. &
            (e2n_local((local_id -1)*numfaces + iface) == -HYDROSTATIC_BC).or.&
            (e2n_local((local_id -1)*numfaces + iface) == -SEEPAGE_BC).or.&
            (e2n_local((local_id -1)*numfaces + iface) == -CONDUCTANCE_BC)) then
          r_p(ghost_face_id) = 0.
        else
          r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface)
        endif
      else
        r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface)
      endif
    enddo
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)

  call PetscLogEventEnd(logging%event_flow_residual_mfd2, ierr);CHKERRQ(ierr)

#endif


end subroutine RichardsResidualPatchMFD2

! ************************************************************************** !

subroutine RichardsResidualPatchMFDLP1(snes,xx,r,realization,ierr)
  ! 
  ! RichardsResidualPatchMFD1: Computes flux and accumulation
  ! terms of the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 
#include "finclude/petscvec.h"
  use petscvec
  
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
  use Material_Aux_class
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), xx_loc_faces_p(:), xx_p(:), work_loc_faces_p(:)

  PetscReal, pointer :: face_fluxes_p(:), bc_faces_p(:)


  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(richards_auxvar_type) :: test_rich_auxvars
  type(global_auxvar_type) :: test_global_auxvars
  class(material_auxvars_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn, bound_id, LP_cell_id 

  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn
  PetscScalar, pointer :: sq_faces(:), e2n_local(:), Smatrix(:,:), face_pr(:)
  PetscScalar, pointer :: neig_den(:), neig_ukvr(:), neig_pres(:), bnd(:), bc_h(:)
  PetscReal :: Res(realization%option%nflowdof), PermTensor(3,3), den, ukvr
  PetscInt :: icell, iface, jface, i,j, numfaces, ghost_face_id, ghost_face_jd
  PetscInt :: ghost_neig_id
  PetscViewer :: viewer

  call PetscLogEventBegin(logging%event_flow_residual_mfd1,  &
                          ierr);CHKERRQ(ierr)
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  discretization => realization%discretization 
  material_auxvars => patch%aux%Material%auxvars
  
  call RichardsUpdateAuxVarsPatchMFDLP(realization)

  patch%aux%Richards%auxvars_up_to_date = PETSC_FALSE ! override flags since they will soon be out of date
  if (option%compute_mass_balance_new) then
!    call RichardsZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr);CHKERRQ(ierr)

  numfaces = 6 ! hex only
  allocate(sq_faces(numfaces))
  allocate(face_pr(numfaces))
  allocate(neig_den(numfaces))
  allocate(neig_ukvr(numfaces))
  allocate(neig_pres(numfaces))
  allocate(bnd(numfaces))
  allocate(bc_h(numfaces))

  do icell = 1, grid%nlmax

    bnd = 0;
    auxvar => grid%MFD%auxvars(icell)

    do j = 1, auxvar%numfaces
      ghosted_id = grid%nL2G(icell)
      ghost_face_id = auxvar%face_id_gh(j)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      jface = grid%faces(ghost_face_id)%id
      sq_faces(j) = conn%area(jface)
      face_pr(j) = xx_loc_faces_p(ghost_face_id)
 
      if ( (e2n_local(j + (icell-1)*numfaces) == -1.0*DIRICHLET_BC).or. &
                (e2n_local(j + (icell-1)*numfaces) == -1.0*HYDROSTATIC_BC))  then
        bnd(j) = 1
      else if ( e2n_local(j + (icell-1)*numfaces) == -1.0*NEUMANN_BC ) then
        bnd(j) = 2
      else if (e2n_local(j + (icell-1)*numfaces)< 0) then
        call printMsg(option,'Type of BC is not supported by MFD mode')
        stop
      endif

      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        neig_pres(j) = xx_loc_faces_p(ghost_face_id)

        if (bnd(j) > 0) then
          call RichardsAuxVarInit(test_rich_auxvars, option)
          call GlobalAuxVarInit(test_global_auxvars, option)

          call RichardsAuxVarCompute(neig_pres(j:j), test_rich_auxvars, &
                       test_global_auxvars, &
                       patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                       material_auxvars(ghosted_id)%porosity, &
                       material_auxvars(ghosted_id)%permeability(perm_xx_index), &
                       option)
        endif

        if (bnd(j)==1) then
          neig_den(j) = test_global_auxvars%den(1)
          neig_ukvr(j) = test_rich_auxvars%kvr
        else if (bnd(j)==2) then
          neig_den(j) = global_auxvars(ghosted_id)%den(1)
        else
          neig_den(j) = global_auxvars(ghosted_id)%den(1)
          neig_ukvr(j) = rich_auxvars(ghosted_id)%kvr
        endif
        if (bnd(j) > 0) call GlobalAuxVarStrip(test_global_auxvars)

      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        if (ghosted_id == conn%id_up(jface)) then
          ghost_neig_id = conn%id_dn(jface)
        else
          ghost_neig_id = conn%id_up(jface)
        endif
 
        neig_den(j) = global_auxvars(ghost_neig_id)%den(1)
        neig_pres(j) = xx_loc_faces_p(grid%ngmax_faces + ghost_neig_id)
#ifdef USE_ANISOTROPIC_MOBILITY
        neig_ukvr(j) = rich_auxvars(ghost_neig_id)%kvr_x
#else
        neig_ukvr(j) = rich_auxvars(ghost_neig_id)%kvr
#endif
      endif
    enddo

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle

    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_id)%permeability(perm_yz_index)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(3,2) = PermTensor(2,3)

    call PetscLogEventBegin(logging%event_flow_flux_mfd, ierr);CHKERRQ(ierr)

    LP_cell_id = grid%ngmax_faces + ghosted_id

    call MFDAuxFluxes(patch, grid, ghosted_id, xx_loc_faces_p(LP_cell_id:LP_cell_id), face_pr, auxvar, PermTensor, &
                      rich_auxvars(ghosted_id), global_auxvars(ghosted_id), &
                      sq_faces,  bnd, neig_den, neig_ukvr, neig_pres, option)
       
    call PetscLogEventEnd(logging%event_flow_flux_mfd, ierr);CHKERRQ(ierr)
      
  enddo

  call VecRestoreArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)

  deallocate(sq_faces)
  deallocate(face_pr)
  deallocate(neig_den)
  deallocate(neig_ukvr)
  deallocate(neig_pres)
  deallocate(bnd)

  call PetscLogEventEnd(logging%event_flow_residual_mfd1, ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsResidualPatchMFDLP1

! ************************************************************************** !

subroutine RichardsResidualPatchMFDLP2(snes,xx,r,realization,ierr)
  ! 
  ! Computes the boundary and source/sink terms of
  ! the residual equation for MFD discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 
#include "finclude/petscvec.h"
  use petscvec
  
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
  use MFD_Aux_module
  use Material_Aux_class
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr

#ifdef DASVYAT

  PetscInt :: i,j
  PetscInt :: local_id, ghosted_id

  PetscScalar, pointer :: r_p(:),  accum_p(:), bc_faces_p(:)
  PetscScalar, pointer :: xx_loc_faces_p(:), flow_xx_p(:)

  PetscScalar :: qsrc, qsrc_mol, tempreal

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:), rich_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(richards_auxvar_type) :: test_rich_auxvars
  type(global_auxvar_type) :: test_global_auxvars
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink, boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection, bc_type, stride


  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn
  PetscScalar :: sq_faces(6), rhs(7), bc_g(6), bc_h(6), face_pres(6), bnd(6)
  PetscScalar :: neig_den(6), neig_kvr(6), neig_dkvr_dp(6), neig_pres(6)
  PetscScalar :: Accum(realization%option%nflowdof), source_f(realization%option%nflowdof)
  PetscScalar :: Res(realization%option%nflowdof), PermTensor(3,3)
  PetscInt :: icell, iface, jface, numfaces, ghost_face_id, ghost_face_jd, ghost_neig_id, LP_cell_id
  PetscScalar, pointer :: e2n_local(:)

  call PetscLogEventBegin(logging%event_flow_residual_mfd2,  &
                          ierr);CHKERRQ(ierr)
  
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

  ! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  numfaces = 6 ! hex only
  stride = 6 !hex only

  do local_id = 1, grid%nlmax

    ghosted_id = grid%nL2G(local_id)
    bc_g = 0
    bc_h = 0
    bnd = 0

    auxvar => grid%MFD%auxvars(local_id)
    do j = 1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(j)
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
        bc_h(j) = bc_faces_p(ghost_face_id)
        bnd(j) = 2
      else if (e2n_local(j + (local_id-1)*stride) < 0) then
        call printMsg(option,'Type of BC is not supported by MFD mode')
        stop
      endif

      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        neig_pres(j) = xx_loc_faces_p(ghost_face_id)

        if (bnd(j) > 0) then
          call RichardsAuxVarInit(test_rich_auxvars, option)
          call GlobalAuxVarInit(test_global_auxvars, option)

          call RichardsAuxVarCompute(neig_pres(j:j), test_rich_auxvars, &
                test_global_auxvars, &
                patch%saturation_function_array(patch%sat_func_id(ghosted_id))%ptr, &
                material_auxvars(ghosted_id)%porosity, &
                materail_auxvars(ghosted_id)%permeaiblity(perm_xx_index), &
                option)
        endif

        if (bnd(j)==1) then
          neig_den(j) = test_global_auxvars%den(1)
          neig_kvr(j) = test_rich_auxvars%kvr
          neig_dkvr_dp(j) =  0.
        else if (bnd(j)==2) then
          neig_den(j) = global_auxvars(ghosted_id)%den(1)
        else
          neig_den(j) = global_auxvars(ghosted_id)%den(1)
          neig_dkvr_dp(j) = rich_auxvars(ghosted_id)%dkvr_dp
          neig_kvr(j) = rich_auxvars(ghosted_id)%kvr
        endif
        if (bnd(j) > 0) call GlobalAuxVarStrip(test_global_auxvars)

      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        if (ghosted_id == conn%id_up(jface)) then
          ghost_neig_id = conn%id_dn(jface)
        else
          ghost_neig_id = conn%id_up(jface)
        endif
 
        neig_den(j) = global_auxvars(ghost_neig_id)%den(1)
        neig_kvr(j) = rich_auxvars(ghost_neig_id)%kvr
        neig_dkvr_dp(j) = rich_auxvars(ghost_neig_id)%dkvr_dp
        neig_pres(j) = xx_loc_faces_p(grid%ngmax_faces + ghost_neig_id)
      endif
    enddo

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsAccumulation(rich_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id)%porosity, &
                              material_auxvars(ghosted_id)%volume, &
                              option, Res)

    Accum(1) = Res(1) - accum_p(local_id)
 
    source_f = 0.
    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_id)%permeability(perm_yz_index)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,2) = PermTensor(2,3)

    call PetscLogEventBegin(logging%event_flow_rhs_mfd, ierr);CHKERRQ(ierr)

    LP_cell_id = grid%ngmax_faces + ghosted_id
         
    call MFDAuxGenerateRhs_LP(patch, grid, ghosted_id, PermTensor, bc_g, source_f, bc_h, auxvar, &
                              rich_auxvars(ghosted_id),&
                              global_auxvars(ghosted_id),&
                              Accum, &
                              material_auxvars(ghosted_id)%porosity, &
                              material_auxvars(ghosted_id)%volume,&
                              xx_loc_faces_p(LP_cell_id:LP_cell_id), face_pres, bnd,&
                              sq_faces, neig_den, neig_kvr, neig_dkvr_dp, neig_pres, option, rhs)

    call PetscLogEventEnd(logging%event_flow_rhs_mfd, ierr);CHKERRQ(ierr)

    do iface=1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(iface)
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
        endif
      else
        r_p(ghost_face_id) = r_p(ghost_face_id) + rhs(iface)
      endif
    enddo

    r_p(grid%ngmax_faces + ghosted_id) = r_p(grid%ngmax_faces + ghosted_id) + rhs(numfaces + 1)

  enddo

  ! Source-sink term
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    if(source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc_mol = qsrc/FMWH2O ! kg/sec -> kmol/sec
        case(HET_MASS_RATE_SS)
          qsrc_mol = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O ! kg/sec -> kmol/sec
        case default
          option%io_buffer='source/sink rate%itype not suppported in MFD'
          call printErrMsg(option)
      end select

      r_p(grid%ngmax_faces+ghosted_id) = r_p(grid%ngmax_faces+ghosted_id) - qsrc_mol
    enddo
    source_sink => source_sink%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_loc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_r_loc_faces, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)

  call PetscLogEventEnd(logging%event_flow_residual_mfd2, ierr);CHKERRQ(ierr)

#endif


end subroutine RichardsResidualPatchMFDLP2

! ************************************************************************** !

subroutine RichardsJacobianMFD(snes,xx,A,B,realization,ierr)
  ! 
  ! RichardsJacobian: Computes the Jacobian for MFD Discretization
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 09/17/10
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
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  
  call PetscLogEventBegin(logging%event_r_jacobian,ierr);CHKERRQ(ierr)

  option => realization%option

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call RichardsJacobianPatchMFD(snes,xx,J,J,realization,ierr)

  if (realization%debug%matview_Jacobian) then
#if 0  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr);CHKERRQ(ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
#endif    
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

  call PetscLogEventEnd(logging%event_r_jacobian,ierr);CHKERRQ(ierr)
  
end subroutine RichardsJacobianMFD

! ************************************************************************** !

subroutine RichardsJacobianMFDLP(snes,xx,A,B,realization,ierr)

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
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  
  call PetscLogEventBegin(logging%event_r_jacobian,ierr);CHKERRQ(ierr)


  option => realization%option

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif


  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call RichardsJacobianPatchMFDLP(snes,xx,J,J,realization,ierr)

  if (realization%debug%matview_Jacobian) then
#if 1
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr);CHKERRQ(ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
#endif    
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

  call PetscLogEventEnd(logging%event_r_jacobian,ierr);CHKERRQ(ierr)

#if 0
    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_dxx_faces.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(realization%field%flow_dxx_faces,viewer,ierr);CHKERRQ(ierr)

    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
 

    call PetscViewerASCIIOpen(realization%option%mycomm,'flow_yy_faces.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(realization%field%flow_yy_faces,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif


end subroutine RichardsJacobianMFDLP

! ************************************************************************** !

subroutine RichardsJacobianPatchMFD (snes,xx,A,B,realization,ierr)
  ! 
  ! RichardsJacobianPatch1: Computes local condensed matrices
  ! for the Jacobian
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 09/17/10
  ! 
       
  
  use MFD_Aux_module
  use Connection_module
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
  use MFD_module
  use Material_Aux_class
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization

  PetscErrorCode :: ierr


#ifdef DASVYAT

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
  PetscReal, pointer :: face_pres(:), xx_faces_p(:), sq_faces(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(richards_parameter_type), pointer :: richards_parameter
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  call VecGetArrayF90(field%flow_xx_loc_faces, xx_faces_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  
  numfaces = 6

  allocate(J(numfaces*numfaces))
  allocate(ghosted_face_id(numfaces))
  allocate(face_pres(numfaces))
  allocate(sq_faces(numfaces))
  allocate(bound_id(numfaces))

  do local_id = 1, grid%nlmax
    auxvar => grid%MFD%auxvars(local_id)

    ghosted_id = grid%nL2G(local_id)
#ifdef USE_ANISOTROPIC_MOBILITY         
    ukvr = rich_auxvars(ghosted_id)%kvr_x
    dkvr_x_dp = rich_auxvars(ghosted_id)%dkvr_x_dp
#else
    ukvr = rich_auxvars(ghosted_id)%kvr
    dkvr_x_dp = rich_auxvars(ghosted_id)%dkvr_dp
#endif
    dden_dp = rich_auxvars(ghosted_id)%dden_dp
    den =  global_auxvars(ghosted_id)%den(1)

    bound_id = 0
    do iface = 1, numfaces
      ghosted_face_id(iface) = auxvar%face_id_gh(iface) - 1
      face_pres(iface) = xx_faces_p(auxvar%face_id_gh(iface))

      conn => grid%faces(ghosted_face_id(iface)+1)%conn_set_ptr
      jface = grid%faces(ghosted_face_id(iface)+1)%id
      sq_faces(iface) = conn%area(jface)

      if ( (e2n_local(iface + (local_id-1)*numfaces) == -DIRICHLET_BC).or. &
            (e2n_local(iface + (local_id-1)*numfaces) == -HYDROSTATIC_BC))  then
        bound_id(iface) = 1
      endif
    enddo

    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_id)%permeability(perm_yz_index)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,2) = PermTensor(2,3)

    J = 0.

    call MFDAuxJacobianLocal( grid, auxvar, &
                              rich_auxvars(ghosted_id), global_auxvars(ghosted_id), &
                              sq_faces, option, J)

    do iface = 1, numfaces
      if (bound_id(iface) == 1) then
        do jface = 1, numfaces
          J((iface-1)*numfaces + jface) = 0.
          J((jface-1)*numfaces + iface) = 0.
        enddo
        J((iface-1)*numfaces + iface) = 1.
      endif
    enddo

    call MatSetValuesLocal(A, numfaces, ghosted_face_id, numfaces, ghosted_face_id, &
                                        J, ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo


  deallocate(ghosted_face_id)
  deallocate(J)
  deallocate(face_pres)
  deallocate(sq_faces)
  deallocate(bound_id)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_mfd.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc_faces, xx_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsJacobianPatchMFD

! ************************************************************************** !

subroutine RichardsJacobianPatchMFDLP (snes,xx,A,B,realization,ierr)
  ! 
  ! RichardsJacobianPatch1: Computes local condensed matrices
  ! for the Jacobian
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 09/17/10
  ! 
       
  
  use MFD_Aux_module
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
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(mfd_auxvar_type), pointer :: auxvar
  type(connection_set_type), pointer :: conn

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  richards_parameter => patch%aux%Richards%richards_parameter
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars

  call VecGetArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)
  numfaces = 6

  allocate(J((numfaces+1)*(numfaces+1)))
  allocate(bound_id(numfaces))
  allocate(sq_faces(numfaces))
  allocate(ghosted_LP_id(numfaces+1+numfaces))
  allocate(neig_LP_id(numfaces))

  ncol = 2*numfaces+1
  nrow = numfaces+1

  do local_id = 1, grid%nlmax
    auxvar => grid%MFD%auxvars(local_id)
    ghosted_id = grid%nL2G(local_id)
    cell_LP_id = grid%ngmax_faces + ghosted_id - 1

    bound_id = 0

    do iface = 1, numfaces
      ghosted_LP_id(iface) = auxvar%face_id_gh(iface) - 1

      conn => grid%faces(ghosted_LP_id(iface)+1)%conn_set_ptr
      jface = grid%faces(ghosted_LP_id(iface)+1)%id
      sq_faces(iface) = conn%area(jface)

      if ( (e2n_local(iface + (local_id-1)*numfaces) == -DIRICHLET_BC).or. &
            (e2n_local(iface + (local_id-1)*numfaces) == -HYDROSTATIC_BC))  then
        bound_id(iface) = 1
      endif
      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        neig_LP_id(iface) = cell_LP_id
      else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
        if (ghosted_id == conn%id_up(jface)) then
          neig_LP_id(iface) = grid%ngmax_faces + conn%id_dn(jface) - 1
        else
          neig_LP_id(iface) = grid%ngmax_faces + conn%id_up(jface) - 1
        endif
      endif
      ghosted_LP_id(numfaces+1+iface) = neig_LP_id(iface)
    enddo
    ghosted_LP_id(numfaces + 1) = grid%ngmax_faces + ghosted_id - 1

    J = 0.

    call MFDAuxJacobianLocal_LP(grid, auxvar, &
                                rich_auxvars(ghosted_id), global_auxvars(ghosted_id), &
                                sq_faces, option, J)
    do jface = 1, numfaces
      if (bound_id(jface) == 1) then
        do iface = 1, numfaces + 1
          J((iface-1)*nrow + jface) = 0.
          J((jface-1)*nrow + iface) = 0.
        enddo
        J(jface + (jface - 1)*nrow) = 1.
      endif
    enddo

    call MatSetValuesLocal(A, numfaces + 1, ghosted_LP_id, numfaces + 1 , ghosted_LP_id, &
                          J, ADD_VALUES,ierr);CHKERRQ(ierr)
    call MatSetValuesLocal(A, 1, cell_LP_id, numfaces, neig_LP_id, auxvar%dRp_dneig, ADD_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  deallocate(J)
  deallocate(bound_id)
  deallocate(ghosted_LP_id)
  deallocate(neig_LP_id)
  deallocate(sq_faces)

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_mfd.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(grid%e2n, e2n_local, ierr);CHKERRQ(ierr)

#endif

end subroutine RichardsJacobianPatchMFDLP

end module Richards_MFD_module
