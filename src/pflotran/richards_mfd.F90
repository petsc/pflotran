module Richards_MFD_module

  use Richards_Aux_module
  use Richards_Common_module
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
  
  public RichardsResidualMFD, RichardsJacobianMFD, &
         RichardsResidualMFDLP, RichardsJacobianMFDLP, &
         RichardsUpdateAuxVarsPatchMFDLP

contains

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
! RichardsResidualPatchMFD1: Computes flux and accumulation 
!   terms of the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsResidualPatchMFD1(snes,xx,r,realization,ierr)

  use Water_EOS_module
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

  use Water_EOS_module
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

  use Water_EOS_module
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

  use Water_EOS_module
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
! RichardsJacobianPatch1: Computes local condensed matrices
!   for the Jacobian
! author: Daniil Svyatskiy
! date: 09/17/10
!
! ************************************************************************** !
subroutine RichardsJacobianPatchMFD (snes,xx,A,B,flag,realization,ierr)
       
  use Water_EOS_module
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
       
  use Water_EOS_module
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

end module Richards_MFD_module
