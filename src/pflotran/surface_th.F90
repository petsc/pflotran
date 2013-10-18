#ifdef SURFACE_FLOW

module Surface_TH_module

  use Surface_Global_Aux_module
  use Surface_TH_Aux_module
  
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
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceTHSetup, &
         SurfaceTHUpdateSurfBC, &
         SurfaceTHUpdateSubsurfSS, &
         SurfaceTHSurf2SubsurfFlux, &
         SurfaceTHCreateSurfSubsurfVec, &
         SurfaceTHGetSubsurfProp, &
         SurfaceTHRHSFunction, &
         SurfaceTHComputeMaxDt, &
         SurfaceTHUpdateAuxVars, &
         SurfaceTHUpdateSolution, &
         SurfaceTHUpdateTemperature

contains

! ************************************************************************** !
!> This routine sets up surface_TH_type
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHSetup(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
 
  implicit none
  
  type(surface_realization_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_aux_vars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_aux_vars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(coupler_type), pointer :: initial_condition
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase
  
  
  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
    
  patch%surf_aux%SurfaceTH => SurfaceTHAuxCreate(option)

  ! allocate aux_var data structures for all grid cells
  allocate(Surf_TH_aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SurfaceTHAuxVarInit(Surf_TH_aux_vars(ghosted_id),option)
  enddo

  patch%surf_aux%SurfaceTH%aux_vars => Surf_TH_aux_vars
  patch%surf_aux%SurfaceTH%num_aux = grid%ngmax

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
    allocate(Surf_TH_aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_aux_vars_bc(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%aux_vars_bc => Surf_TH_aux_vars_bc
  endif
  patch%surf_aux%SurfaceTH%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (sum_connection > 0) then
    allocate(Surf_TH_aux_vars_ss(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_aux_vars_ss(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%aux_vars_ss => Surf_TH_aux_vars_ss
  endif
  patch%surf_aux%SurfaceTH%num_aux_ss = sum_connection

  call SurfaceTHSetPlotVariables(surf_realization)

end subroutine SurfaceTHSetup

! ************************************************************************** !
!> This routine adds variables to be printed to list
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHSetPlotVariables(surf_realization)
  
  use Surface_Realization_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => surf_realization%output_option%output_variable_list
  
  name = 'H'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_HEAD)

  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_TEMPERATURE)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceTHSetPlotVariables

! ************************************************************************** !
!> This routine gets latest states (P,T) from subsurface model and updates
!! boundary condition for surface flow model.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHUpdateSurfBC(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module
  use Global_Aux_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  type(global_auxvar_type), pointer :: global_aux_vars(:)  
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: qsrc_p(:),press_p(:),temp_p(:)
  PetscInt :: local_id,iconn,sum_connection,ghosted_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: yy_p(:)
  PetscReal, pointer :: zz_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: dist_p(:)
  PetscReal :: dist_x,dist_y,dist_z,dist
  PetscInt :: size,istart,iend
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field
  global_aux_vars => patch%aux%Global%aux_vars

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
  
  ! Update the surface BC
  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  sum_connection = 0
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      cur_connection_set => coupler%connection_set
      
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then

        ! Exchange subsurface PRESSURE
        call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,press_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          iend = ghosted_id*option%nflowdof
          istart = iend-option%nflowdof+1
          press_p(iconn) = xx_loc_p(istart)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,press_p,ierr)
        call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)

        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! Exchange subsurface TEMPERATURE
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,temp_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          iend = ghosted_id*option%nflowdof
          istart = iend-option%nflowdof+1
          temp_p(iconn) = global_aux_vars(ghosted_id)%temp(1)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,temp_p,ierr)

        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%temp_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                        surf_field%subsurf_temp_vec_1dof,surf_field%temp_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)

      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    coupler => coupler%next
  enddo

end subroutine SurfaceTHUpdateSurfBC

! RTM: TODO: Figure out if this needs to be modified for surface freezing.
! ************************************************************************** !
!> This routine updates source/sink term for the subsurface model
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHUpdateSubsurfSS(realization,surf_realization,dt)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_loc_p(:),vec_p(:)
  PetscReal :: dt
  PetscInt :: local_id
  PetscInt :: iconn
  PetscReal :: den
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,NFLOWDOF)

  call density(option%reference_temperature,option%reference_pressure,den)

  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then
      
        coupler_found = PETSC_TRUE
        
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_ndof, &
                             surf_field%exchange_subsurf_2_surf, &
                             surf_field%subsurf_temp_vec_ndof, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_ndof, &
                           surf_field%exchange_subsurf_2_surf, &
                           surf_field%subsurf_temp_vec_ndof, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        call VecGetArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)
        do iconn=1,coupler%connection_set%num_connections
          ! Flux of water
          coupler%flow_aux_real_var(ONE_INTEGER,iconn)=-vec_p((iconn-1)*option%nflowdof+1)/dt*den
          
          ! Heat flux
          coupler%flow_aux_real_var(TWO_INTEGER,iconn)=-vec_p(iconn*option%nflowdof)/dt
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)

        call VecSet(surf_field%exchange_subsurf_2_surf,0.d0,ierr)
      endif

    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_ss.'
    call printErrMsg(option)
  endif
  
end subroutine SurfaceTHUpdateSubsurfSS

! ************************************************************************** !
!> This routine creates a PETSc vector to data between surface and subsurface
!! model.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHCreateSurfSubsurfVec(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_loc_p(:),vec_p(:)
  PetscReal :: den          ! density      [kg/m^3]
  PetscInt :: local_id,i
  PetscInt :: iconn
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then
        if(coupler%flow_condition%rate%itype/=HET_MASS_RATE_SS) then
          option%io_buffer = 'Flow condition from_surface_ss should be ' // &
            'heterogeneous_mass_rate'
          call printErrMsg(option)
        endif
        coupler_found = PETSC_TRUE

        call VecCreate(option%mycomm,surf_field%subsurf_temp_vec_ndof,ierr)
        call VecSetSizes(surf_field%subsurf_temp_vec_ndof, &
                         coupler%connection_set%num_connections*option%nflowdof, &
                         PETSC_DECIDE,ierr)
        call VecSetFromOptions(surf_field%subsurf_temp_vec_ndof,ierr)

        call VecCreate(option%mycomm,surf_field%subsurf_temp_vec_1dof,ierr)
        call VecSetSizes(surf_field%subsurf_temp_vec_1dof, &
                         coupler%connection_set%num_connections, &
                         PETSC_DECIDE,ierr)
        call VecSetFromOptions(surf_field%subsurf_temp_vec_1dof,ierr)

      endif
    endif
    
    coupler => coupler%next
  enddo

  if(.not.coupler_found) then
    option%io_buffer = 'Missing within the input deck for subsurface ' // &
      'boundary condition named from_surface_ss.'
    call printErrMsg(option)
  endif

end subroutine SurfaceTHCreateSurfSubsurfVec

! ************************************************************************** !
!> This routine computes the flux of water and energy from surface to 
!! subsurface model.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHSurf2SubsurfFlux(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Water_EOS_module
  use Saturation_Function_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module
  use Surface_TH_Aux_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars(:)
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal :: den_surf_kg          ! density      [kg/m^3]
  PetscReal, parameter :: den_surf_ice_kg=917.d0   ! density of surface ice [kg/m^3]
  PetscReal :: den_sub_kg           ! density      [kg/m^3]
  PetscReal :: den_aveg             ! density      [kg/m^3]
  PetscInt :: local_id, ghosted_id, iconn
  
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: press_sub_p(:) ! Pressure [Pa]
  PetscReal, pointer :: temp_sub_p(:)  ! Temperature [C]
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: exch_p(:)
  PetscReal, pointer :: area_p(:)
  PetscReal, pointer :: dist_gravity_p(:)
  PetscReal, pointer :: dist_p(:)
  PetscReal, pointer :: ckdry_p(:)
  PetscReal, pointer :: ckwet_p(:)
  PetscReal, pointer :: ckice_p(:)
  PetscReal, pointer :: th_alpha_p(:)
  PetscReal, pointer :: th_alpha_fr_p(:)
  PetscReal, pointer :: sat_ice_p(:)
  
  PetscReal :: hw
  PetscReal :: press_surf
  PetscReal :: dphi
  PetscReal :: press
  PetscReal :: sat
  PetscReal :: kr
  PetscReal :: ds_dp
  PetscReal :: dkr_dp
  PetscReal :: por
  PetscReal :: perm
  PetscBool :: saturated
  PetscReal :: sat_pressure
  PetscReal :: pw
  PetscReal :: visl
  PetscReal :: dvis_dp
  PetscReal :: dvis_dt
  PetscReal :: v_darcy
  PetscReal :: v_darcy_max
  PetscReal :: gravity
  PetscReal :: press_up, press_dn
  PetscReal :: k_eff_dn, k_eff_up
  PetscReal :: Dk_eff
  PetscReal :: Ke_up, Ke_fr, Ke_fr_up
  PetscReal :: dtemp
  PetscReal :: Cwi
  PetscReal :: temp_half
  PetscReal, parameter :: epsilon = 1.d-6
    
  PetscBool :: coupler_found = PETSC_FALSE
  PetscBool :: v_darcy_limit

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field
  surf_global_aux_vars => surf_patch%surf_aux%SurfaceGlobal%aux_vars
  surf_aux_vars => surf_patch%surf_aux%SurfaceTH%aux_vars

  call VecGetArrayF90(surf_field%press_subsurf,press_sub_p,ierr)
  call VecGetArrayF90(surf_field%temp_subsurf,temp_sub_p,ierr)
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_p,ierr)
  call VecGetArrayF90(surf_field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(surf_field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(surf_field%Dq,Dq_p,ierr)
  call VecGetArrayF90(surf_field%exchange_subsurf_2_surf,exch_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist_gravity,dist_gravity_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist,dist_p,ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)
  call VecGetArrayF90(surf_field%ckdry,ckdry_p,ierr)
  call VecGetArrayF90(surf_field%ckwet,ckwet_p,ierr)
  call VecGetArrayF90(surf_field%ckice,ckice_p,ierr)
  call VecGetArrayF90(surf_field%th_alpha,th_alpha_p,ierr)
  call VecGetArrayF90(surf_field%th_alpha_fr,th_alpha_fr_p,ierr)
  call VecGetArrayF90(surf_field%sat_ice,sat_ice_p,ierr)

  ! Update the surface BC
  coupler_list => surf_patch%source_sinks
  coupler => coupler_list%first

  v_darcy_max=0.d0
  v_darcy_limit=PETSC_FALSE
  
  !
  !            SURFACE
  !              (DN)
  !   ---------------------------------
  !           SUBSURFACE         ///\\\
  !              (UP)
  !

  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if(StringCompare(coupler%name,'from_subsurface_ss')) then

      if (.not.associated(coupler%flow_aux_real_var)) then
        option%io_buffer='SURF_SOURCE_SINK: from_subsurface_ss was does not ' //&
          ' a heterogeneous flow condition associated with it'
        call printErrMsg(option)
      endif

      if(coupler%flow_condition%rate%itype/=HET_VOL_RATE_SS) then
        option%io_buffer = 'Flow condition from_subsurface_ss should be ' // &
          'heterogeneous_volumetric_rate'
        call printErrMsg(option)
      endif
      
      cur_connection_set => coupler%connection_set

      do iconn = 1, cur_connection_set%num_connections

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = surf_grid%nL2G(local_id)

        ! Compute densities:
        call density(surf_global_aux_vars(ghosted_id)%temp(1), &
                     option%reference_pressure,den_surf_kg)
        ! Now modify den_surf_kg to account for frozen fraction.
        ! WARNING: This assumes density of ice at atmospheric pressure;
        ! TODO: Need to actually compute this to handle the general case.
        den_surf_kg = surf_aux_vars(ghosted_id)%unfrozen_fraction*den_surf_kg + &
                      (1-surf_aux_vars(ghosted_id)%unfrozen_fraction)*den_surf_ice_kg
        call density(temp_sub_p(local_id),press_sub_p(local_id),den_sub_kg)
        den_aveg = (den_surf_kg + den_sub_kg)/2.d0

        ! Exchange of water between surface-subsurface.
        hw = xx_p((local_id-1)*option%nflowdof+1)
        press_surf = hw*(abs(option%gravity(3)))*den_surf_kg + &
                     option%reference_pressure

        press_up = press_sub_p(local_id)
        press_dn = press_surf
        gravity = dist_gravity_p(local_id)*den_aveg
        
        dphi = press_up - press_dn + gravity
        
        ! check if there is standing water on the surface
        if (dphi<0.d0 .and. press_dn - option%reference_pressure<eps) then
          dphi= 0.d0
        endif
        
        ! exfiltration will only occur if subsurface pressure is greater
        ! than reference pressure
        if (dphi>=0 .and. press_up-option%reference_pressure<eps ) then
          dphi = 0.d0
        endif
        
        if(dphi>0.d0) then
          press = press_up
        else
          press = press_dn
        endif
        
        call SaturationFunctionCompute( &
          press,sat,kr,ds_dp,dkr_dp,&
          patch%saturation_function_array(int(icap_loc_p(local_id)))%ptr, &
          0.d0,0.d0,saturated,option)
        
        if(saturated) then
          pw = press
        else
          pw = option%reference_pressure
        endif
                                           
        call psat(option%reference_temperature,sat_pressure,ierr)
        call VISW(option%reference_temperature,pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr)

        v_darcy = Dq_p(local_id)*kr/visl*dphi
        if (v_darcy<=0.d0) then
          ! Flow is happening from surface to subsurface.
          ! Note that we limit the mass exchange when surface ice is present 
          ! (as the frozen fraction is immobile) by multiplying by the 
          ! unfrozen fraction.
          if ( abs(v_darcy) > surf_aux_vars(ghosted_id)%unfrozen_fraction * &
                              xx_p(local_id)/option%surf_flow_dt ) then
            v_darcy = surf_aux_vars(ghosted_id)%unfrozen_fraction * &
                      (-xx_p(local_id)/option%surf_flow_dt)
            v_darcy_limit=PETSC_TRUE
          endif
          temp_half = surf_global_aux_vars(ghosted_id)%temp(1) + 273.15d0
        else
          ! Exfiltration is occuring
          temp_half = temp_sub_p(local_id) + 273.15d0
        endif
        
        ! Mass flux
        exch_p((local_id-1)*option%nflowdof+1) = &
          exch_p((local_id-1)*option%nflowdof+1) + &
          v_darcy*area_p(local_id)*option%surf_flow_dt
        !coupler%flow_aux_real_var(ONE_INTEGER,local_id)=v_darcy
        coupler%flow_aux_real_var(ONE_INTEGER,local_id)=0.d0
        xx_p((local_id-1)*option%nflowdof+1) = &
          xx_p((local_id-1)*option%nflowdof+1) + v_darcy*option%surf_flow_dt
        if(abs(v_darcy)>v_darcy_max) v_darcy_max = v_darcy

        ! Heat flux associated with mass flux
        Cwi = surf_aux_vars(ghosted_id)%Cwi
        exch_p((local_id-1)*option%nflowdof+2) = &
          exch_p((local_id-1)*option%nflowdof+2) + &
          den_aveg*v_darcy*temp_half*Cwi*area_p(local_id)*option%surf_flow_dt

        if (hw>0.d0) then
          ! Exchange of heat between surface water--subsurface domain
          Ke_up = (sat + epsilon)**th_alpha_p(local_id)
          k_eff_up = ckdry_p(local_id) + &
                      (ckwet_p(local_id) - ckdry_p(local_id))*Ke_up
#ifdef ICE
          Ke_fr_up = (sat_ice_p(local_id) + epsilon)**th_alpha_fr_p(local_id)
          k_eff_up = ckwet_p(local_id)*Ke_up + ckice_p(local_id)*Ke_fr_up + &
                     ckdry_p(local_id)*(1.d0 - Ke_up - Ke_fr_up)
#endif
          k_eff_dn = surf_aux_vars(ghosted_id)%k_therm

          Dk_eff = k_eff_up*k_eff_dn/(k_eff_up*hw/2.d0 + &
                                      k_eff_dn*dist_p(local_id))

          dtemp = temp_sub_p(local_id) - surf_global_aux_vars(ghosted_id)%temp(1)
          exch_p(local_id*option%nflowdof) = exch_p(local_id*option%nflowdof) + &
            area_p(local_id)*Dk_eff*dTemp*option%surf_flow_dt
        endif
      enddo
    else

      ! In the absence of standing water, the heat flux SS term for surface flow
      ! is directly applied to subsurface domain.

      if (associated(coupler%flow_aux_real_var) .and. &
          coupler%flow_condition%energy_rate%itype == HET_ENERGY_RATE_SS) then

        cur_connection_set => coupler%connection_set

        do iconn = 1, cur_connection_set%num_connections

          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = surf_grid%nL2G(local_id)

          hw = xx_p((local_id-1)*option%nflowdof+1)

          if (hw==0.d0) then
            exch_p(local_id*option%nflowdof) = exch_p(local_id*option%nflowdof) &
              -coupler%flow_aux_real_var(TWO_INTEGER,local_id)*option%surf_flow_dt
            ! NOTE: There is a -ve sign
            !   coupler%flow_aux_real_var(TWO_INTEGER,:) : +ve value implies a
            !       source of energy to surface water; while
            !   exch_p(:) : +ve values implies transfer of energy from subsurface
            !       to surface water.
          endif

        enddo

      endif

    endif
    coupler => coupler%next
  enddo
  
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)
  call VecRestoreArrayF90(surf_field%surf2subsurf_dist_gravity,dist_gravity_p,ierr)
  call VecRestoreArrayF90(surf_field%surf2subsurf_dist,dist_p,ierr)
  call VecRestoreArrayF90(surf_field%exchange_subsurf_2_surf,exch_p,ierr)
  call VecRestoreArrayF90(surf_field%Dq,Dq_p,ierr)  
  call VecRestoreArrayF90(surf_field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_p,ierr)
  call VecRestoreArrayF90(surf_field%temp_subsurf,temp_sub_p,ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf,press_sub_p,ierr)
  call VecRestoreArrayF90(surf_field%ckdry,ckdry_p,ierr)
  call VecRestoreArrayF90(surf_field%ckwet,ckwet_p,ierr)
  call VecRestoreArrayF90(surf_field%ckice,ckice_p,ierr)
  call VecRestoreArrayF90(surf_field%th_alpha,th_alpha_p,ierr)
  call VecRestoreArrayF90(surf_field%th_alpha_fr,th_alpha_fr_p,ierr)
  call VecRestoreArrayF90(surf_field%sat_ice,sat_ice_p,ierr)

end subroutine SurfaceTHSurf2SubsurfFlux

! ************************************************************************** !
!> This routine get soil properties of the top-most soil layer from the 
!! subsurface domain.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 02/28/13
! ************************************************************************** !
subroutine SurfaceTHGetSubsurfProp(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Field_module
  use Water_EOS_module
  use Discretization_module
  use Connection_module
  use Surface_Realization_class
  use Realization_Base_class
  use DM_Kludge_module
  use TH_Aux_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type)         :: realization
  type(surface_realization_type) :: surf_realization

  type(patch_type),pointer            :: patch,surf_patch
  type(grid_type),pointer             :: grid,surf_grid
  type(coupler_list_type), pointer    :: coupler_list
  type(coupler_type), pointer         :: coupler
  type(flow_condition_type), pointer  :: flow_condition
  type(option_type), pointer          :: option
  type(field_type),pointer            :: field
  type(surface_field_type),pointer    :: surf_field
  type(dm_ptr_type), pointer          :: dm_ptr
  type(connection_set_type), pointer  :: cur_connection_set
  type(TH_parameter_type), pointer :: TH_parameter
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal, pointer :: qsrc_p(:),vec_p(:)
  PetscInt :: local_id,iconn,sum_connection,ghosted_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: yy_p(:)
  PetscReal, pointer :: zz_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: dist_gravity_p(:)
  PetscReal, pointer :: dist_p(:)
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: ckwet_p(:)
  PetscReal, pointer :: ckdry_p(:)
  PetscReal, pointer :: ckice_p(:)
  PetscReal, pointer :: th_alpha_p(:)
  PetscReal, pointer :: th_alpha_fr_p(:)
  PetscReal :: dist_x,dist_y,dist_z,dist
  PetscInt :: size
  
  PetscBool :: coupler_found = PETSC_FALSE

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field
  TH_parameter => patch%aux%TH%TH_parameter

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
  
  ! Update the surface BC
  coupler_list => patch%source_sinks
  coupler => coupler_list%first
  sum_connection = 0
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then
      cur_connection_set => coupler%connection_set
      
      ! Find the BC from the list of BCs
      if(StringCompare(coupler%name,'from_surface_ss')) then

        ! perm_x
        call VecGetArrayF90(field%perm_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_xx_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                             surf_field%subsurf_temp_vec_1dof, &
                             surf_field%perm_xx, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_xx, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! perm_y
        call VecGetArrayF90(field%perm_yy_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_yy_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%perm_yy, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_yy, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! perm_z
        call VecGetArrayF90(field%perm_zz_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%perm_zz_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%perm_zz, &
                            INSERT_VALUES, SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%perm_zz, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! por
        call VecGetArrayF90(field%porosity_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(field%porosity_loc,xx_loc_p, ierr)
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%por, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%por, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! icap
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=icap_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%icap_loc, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%icap_loc, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! ithrm_id
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=ithrm_loc_p(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%ithrm_loc, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%ithrm_loc, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! x
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%x(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_xx, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_xx, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! y
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%y(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_yy, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_yy, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        ! z
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=grid%z(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%subsurf_zz, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%subsurf_zz, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    coupler => coupler%next
  enddo

  call VecGetArrayF90(surf_field%subsurf_xx,xx_p,ierr)
  call VecGetArrayF90(surf_field%subsurf_yy,yy_p,ierr)
  call VecGetArrayF90(surf_field%subsurf_zz,zz_p,ierr)
  call VecGetArrayF90(surf_field%perm_xx,perm_xx_p,ierr)
  call VecGetArrayF90(surf_field%perm_yy,perm_yy_p,ierr)
  call VecGetArrayF90(surf_field%perm_zz,perm_zz_p,ierr)
  call VecGetArrayF90(surf_field%Dq,Dq_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist_gravity,dist_gravity_p,ierr)
  call VecGetArrayF90(surf_field%surf2subsurf_dist,dist_p,ierr)
  call VecGetArrayF90(surf_field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(surf_field%ckwet,ckwet_p,ierr)
  call VecGetArrayF90(surf_field%ckdry,ckdry_p,ierr)
  call VecGetArrayF90(surf_field%ckice,ckice_p,ierr)
  call VecGetArrayF90(surf_field%th_alpha,th_alpha_p,ierr)
  call VecGetArrayF90(surf_field%th_alpha_fr,th_alpha_fr_p,ierr)

  do local_id=1,surf_grid%nlmax
    dist_x = (xx_p(local_id) - surf_grid%x(local_id))
    dist_y = (yy_p(local_id) - surf_grid%y(local_id))
    dist_z = (zz_p(local_id) - surf_grid%z(local_id))
      
    dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z)

    Dq_p(local_id) = (perm_xx_p(local_id)*abs(dist_x)/dist + &
                      perm_yy_p(local_id)*abs(dist_y)/dist + &
                      perm_zz_p(local_id)*abs(dist_z)/dist)/dist
    dist_gravity_p(local_id) = dist*(dist_x*option%gravity(1)+ &
                             dist_y*option%gravity(2)+ &
                             dist_z*option%gravity(3))

    dist_p(local_id) = dist
    ckwet_p(local_id) = TH_parameter%ckwet(int(ithrm_loc_p(local_id)))
    ckdry_p(local_id) = TH_parameter%ckdry(int(ithrm_loc_p(local_id)))
    th_alpha_p(local_id) = TH_parameter%alpha(int(ithrm_loc_p(local_id)))

#ifdef ICE
    ckice_p(local_id) = TH_parameter%ckfrozen(int(ithrm_loc_p(local_id)))
    th_alpha_fr_p(local_id) = TH_parameter%alpha_fr(int(ithrm_loc_p(local_id)))
#else
    ckice_p(local_id) = 0.d0
    th_alpha_fr_p(local_id) = 0.d0
#endif

  enddo

  call VecRestoreArrayF90(surf_field%surf2subsurf_dist_gravity,dist_gravity_p,ierr)
  call VecRestoreArrayF90(surf_field%surf2subsurf_dist,dist_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_xx,xx_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_yy,yy_p,ierr)
  call VecRestoreArrayF90(surf_field%subsurf_zz,zz_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_xx,perm_xx_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_yy,perm_yy_p,ierr)
  call VecRestoreArrayF90(surf_field%perm_zz,perm_zz_p,ierr)
  call VecRestoreArrayF90(surf_field%Dq,Dq_p,ierr)
  call VecRestoreArrayF90(surf_field%ckwet,ckwet_p,ierr)
  call VecRestoreArrayF90(surf_field%ckdry,ckdry_p,ierr)
  call VecRestoreArrayF90(surf_field%ckice,ckice_p,ierr)
  call VecRestoreArrayF90(surf_field%th_alpha,th_alpha_p,ierr)
  call VecRestoreArrayF90(surf_field%th_alpha_fr,th_alpha_fr_p,ierr)

  surf_realization%first_time=PETSC_FALSE

end subroutine SurfaceTHGetSubsurfProp

! ************************************************************************** !
!> This routine provides the function evaluation for PETSc TSSolve()
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceTHRHSFunction(ts,t,xx,ff,surf_realization,ierr)

  use Water_EOS_module
  use Connection_module
  use Surface_Realization_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  TS                             :: ts
  PetscReal                      :: t
  Vec                            :: xx
  Vec                            :: ff
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(coupler_type), pointer               :: source_sink
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set

  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_ss(:)

  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt :: istart, iend

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: qsrc, qsrc_flow
  PetscReal :: esrc

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH)       :: string,string2

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  surf_aux_vars => patch%surf_aux%SurfaceTH%aux_vars
  surf_aux_vars_bc => patch%surf_aux%SurfaceTH%aux_vars_bc
  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc
  surf_global_aux_vars_ss => patch%surf_aux%SurfaceGlobal%aux_vars_ss

  write(string,*),t
  call printMsg(option,"In SurfaceTHRHSFunction: " // string)
  
  surf_realization%iter_count = surf_realization%iter_count+1
  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Then, update the aux vars
  ! RTM: This includes calculation of the accumulation terms, correct?
  call SurfaceTHUpdateAuxVars(surf_realization)
  ! override flags since they will soon be out of date  
  patch%surf_aux%SurfaceTH%aux_vars_up_to_date = PETSC_FALSE

  call VecGetArrayF90(ff,ff_p, ierr)
  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p, ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  ff_p = 0.d0
  Res  = 0.d0

  ! RTM: Does this computed density get used anywhere?
  call density(option%reference_temperature,option%reference_pressure,rho)

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

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

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
#if 0
      call SurfaceTHFlux(surf_aux_vars(ghosted_id_up), &
                         surf_global_aux_vars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_aux_vars(ghosted_id_dn), &
                         surf_global_aux_vars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,Res)
#endif

      !patch%internal_velocities(1,sum_connection) = vel
      !patch%surf_internal_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(TH_PRESSURE_DOF)
      !patch%surf_internal_fluxes(TH_TEMPERATURE_DOF,sum_connection) =
      vel = patch%internal_velocities(1,sum_connection)
      Res(:) = patch%surf_internal_fluxes(:,sum_connection)

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) - Res(:)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
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
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

#if 0
      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         surf_aux_vars_bc(sum_connection), &
                         surf_global_aux_vars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)
#endif

      !patch%boundary_velocities(1,sum_connection) = vel
      !patch%surf_boundary_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(1)

      vel = patch%boundary_velocities(1,sum_connection)
      Res(:) = patch%surf_boundary_fluxes(:,sum_connection)
      
      iend = local_id_dn*option%nflowdof
      istart = iend-option%nflowdof+1
      ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    if(source_sink%flow_condition%rate%itype/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%dataset%rarray(1)
      
    if(source_sink%flow_condition%rate%itype/=HET_ENERGY_RATE_SS) &
      esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc = m^3/sec
          qsrc = qsrc_flow*area_p(local_id)
        case(HET_VOL_RATE_SS)
          ! qsrc = m^3/sec
          qsrc = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)*area_p(local_id)
        case default
          option%io_buffer = 'Source/Sink flow condition type not recognized'
          call printErrMsg(option)
      end select
      
      esrc = 0.d0
      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (ENERGY_RATE_SS)
          esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (HET_ENERGY_RATE_SS)
          esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ff_p(istart) = ff_p(istart) + qsrc/area_p(local_id)
      ! RTM: TODO: What should the density term and specific heat capactiy be 
      ! in the freezing case?
      ! I think using the weighted average of liquid and ice densities and Cwi 
      ! is correct here, but I should check.
      ff_p(iend) = ff_p(iend) + esrc + &
                    surf_global_aux_vars_ss(local_id)%den_kg(1)* &
                    (surf_global_aux_vars_ss(local_id)%temp(1) + 237.15d0)* &
                    surf_aux_vars(local_id)%Cwi* &
                    qsrc/area_p(local_id)* &
                    qsrc*option%surf_flow_dt

    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(ff,ff_p, ierr)
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)

  if (surf_realization%debug%vecview_solution) then
    string = 'Surf_xx_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    string = 'Surf_ff_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr)
    call VecView(ff,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine SurfaceTHRHSFunction

! ************************************************************************** !
!> This routine maximum allowable 'dt' for explicit time scheme.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date:
! ************************************************************************** !
subroutine SurfaceTHComputeMaxDt(surf_realization,max_allowable_dt)

  use Water_EOS_module
  use Connection_module
  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  PetscErrorCode                 :: ierr

  type(grid_type), pointer                  :: grid
  type(patch_type), pointer                 :: patch
  type(option_type), pointer                :: option
  type(surface_field_type), pointer         :: surf_field
  type(coupler_type), pointer               :: boundary_condition
  type(connection_set_list_type), pointer   :: connection_set_list
  type(connection_set_type), pointer        :: cur_connection_set

  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)

  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: max_allowable_dt

  PetscReal, pointer :: mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  surf_aux_vars => patch%surf_aux%SurfaceTH%aux_vars
  surf_aux_vars_bc => patch%surf_aux%SurfaceTH%aux_vars_bc
  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc

  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p, ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10
  vel = 0.d0

  call density(option%reference_temperature,option%reference_pressure,rho)
  !call nacl_den(option%reference_temperature,option%reference_pressure*1d-6,0.d0,rho)
  !rho = rho * 1.d3

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

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

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      call SurfaceTHFlux(surf_aux_vars(ghosted_id_up), &
                         surf_global_aux_vars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_aux_vars(ghosted_id_dn), &
                         surf_global_aux_vars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%internal_velocities(1,sum_connection) = vel
      patch%surf_internal_fluxes(:,sum_connection) = Res(:)
      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/4.d0)

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
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         surf_aux_vars_bc(sum_connection), &
                         surf_global_aux_vars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%surf_boundary_fluxes(:,sum_connection) = Res(:)

      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/4.d0)
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p,ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr)
  
end subroutine SurfaceTHComputeMaxDt

! ************************************************************************** !
!> This routine computes the internal flux term for under
!! diffusion-wave assumption.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceTHFlux(surf_aux_var_up, &
                         surf_global_aux_var_up, &
                         zc_up, &
                         mannings_up, &
                         surf_aux_var_dn, &
                         surf_global_aux_var_dn, &
                         zc_dn, &
                         mannings_dn, &
                         dist, &
                         length, &
                         option, &
                         vel, &
                         Res)

  use Surface_TH_Aux_module
  use Surface_Global_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_aux_var_up
  type(Surface_TH_auxvar_type) :: surf_aux_var_dn
  type(surface_global_auxvar_type) :: surf_global_aux_var_up
  type(surface_global_auxvar_type) :: surf_global_aux_var_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn

  PetscReal :: head_up, head_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s
  
  PetscReal :: flux       ! units: m^2/s
  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: unfrozen_fraction_half
  PetscReal :: dhead
  PetscReal :: den_aveg
  PetscReal :: temp_half
  PetscReal :: dtemp
  PetscReal :: Cw
  PetscReal :: k_therm

  ! initialize
  flux = 0.d0

  ! Flow equation
  head_up = surf_global_aux_var_up%head(1) + zc_up
  head_dn = surf_global_aux_var_dn%head(1) + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    temp_half = surf_global_aux_var_up%temp(1) + 273.15d0
    unfrozen_fraction_half = surf_aux_var_up%unfrozen_fraction
    if (surf_global_aux_var_up%head(1)>eps) then
      hw_half = surf_global_aux_var_up%head(1)
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    temp_half = surf_global_aux_var_dn%temp(1) + 273.15d0
    unfrozen_fraction_half = surf_aux_var_dn%unfrozen_fraction
    if (surf_global_aux_var_dn%head(1)>eps) then
      hw_half = surf_global_aux_var_dn%head(1)
    else
      hw_half = 0.d0
    endif
  endif
  
  dhead=head_up-head_dn
  if(abs(dhead)<eps) then
    dhead=0.d0
    vel = 0.d0
  else
    ! RTM: We modify the term raised to the power 2/3 (the "hydraulic radius") 
    ! by the (upwinded) unfrozen fraction.  For a wide rectangular channel, 
    ! hydraulic radius (which is a measure of the "efficiency" of the channel) 
    ! is often taken to be the flow depth, so I believe this makes sense. (?)
    ! The actual total head term ('hw_half' here) is NOT modified by the 
    ! unfrozen fraction: though the ice is immobile, its weight does 
    ! contribute to the pressure head.
    vel = ((unfrozen_fraction_half * hw_half)**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)

     !RTM: Original code for when freezing is not considered is 
!    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
!          dhead/(abs(dhead)**(1.d0/2.d0))* &
!          1.d0/(dist**0.5d0)
  endif

  flux = unfrozen_fraction_half * hw_half*vel
  Res(TH_PRESSURE_DOF) = flux*length
  
  ! Temperature equation
  ! RTM: k_therm is the weighted average of the liquid and ice thermal 
  ! conductivities.  For the density and specific heat capacity in the 
  ! advection term, we want these for liquid water ONLY, as the ice portion 
  ! is immobile and thus should not make up part of the advection term. We 
  ! also multiply the ponded water depth (hw_half) by the unfrozen fraction 
  ! in the advection term but NOT the conduction term.
  ! We do the same in SurfaceTHBCFlux().

  ! Average density
  ! Here we only consider the LIQUID fraction.
  den_aveg = (surf_aux_var_up%den_water_kg + &
              surf_aux_var_dn%den_water_kg)/2.d0
  ! Temperature difference
  dtemp = surf_global_aux_var_up%temp(1) - surf_global_aux_var_dn%temp(1)

  ! Note, Cw and k_therm are same for up and downwind
  Cw = surf_aux_var_up%Cw
  k_therm = surf_aux_var_up%k_therm
  
  ! Unfrozen fraction multiplies hw_half in advection term, but does NOT affect the 
  ! conduction therm.  
  ! RTM: Brookfield et al. 2009 also has dispersion term, which we are not using.
  Res(TH_TEMPERATURE_DOF) = (den_aveg*vel*temp_half*Cw*unfrozen_fraction_half*hw_half + &
                             k_therm*dtemp/dist*hw_half)*length

end subroutine SurfaceTHFlux

! ************************************************************************** !
!> This routine computes flux for boundary cells.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHBCFlux(ibndtype, &
                           surf_aux_var, &
                           surf_global_aux_var, &
                           slope, &
                           mannings, &
                           length, &
                           option, &
                           vel, &
                           Res)

  use Option_module
  
  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_aux_var
  type(surface_global_auxvar_type) :: surf_global_aux_var
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt  :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: Res(1:option%nflowdof) 

  PetscInt :: pressure_bc_type
  PetscReal :: head

  flux = 0.d0
  vel = 0.d0
  
  ! RTM: I've multiplied the head (ponded water depth, actually) by the 
  ! unfrozen fraction.  I believe this makes sense, but I should think a bit 
  ! more about what a "zero gradient" condition means in the case of freezing
  ! surface water.

  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_aux_var%head(1)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
      else
        vel = -sqrt(dabs(slope))/mannings*((surf_aux_var%unfrozen_fraction * head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = surf_aux_var%unfrozen_fraction * head*vel
  Res(TH_PRESSURE_DOF) = flux*length

  ! Temperature
  ! RTM: See note about in SufaceTHFlux() about how frozen/unfrozen are handled here.
  Res(TH_TEMPERATURE_DOF) = surf_aux_var%den_water_kg* &
                            (surf_global_aux_var%temp(1) + 273.15d0)* &
                            surf_aux_var%Cw* &
                            vel*head*surf_aux_var%unfrozen_fraction*length

end subroutine SurfaceTHBCFlux

! ************************************************************************** !
!> This routine updates auxiliary variables
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHUpdateAuxVars(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_th_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_aux_vars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_aux_vars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_th_aux_vars => patch%surf_aux%SurfaceTH%aux_vars
  surf_th_aux_vars_bc => patch%surf_aux%SurfaceTH%aux_vars_bc
  surf_th_aux_vars_ss => patch%surf_aux%SurfaceTH%aux_vars_ss
  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc
  surf_global_aux_vars_ss => patch%surf_aux%SurfaceGlobal%aux_vars_ss
  
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call SurfaceTHAuxVarCompute(xx_loc_p(istart:iend), &
                                surf_th_aux_vars(ghosted_id), &
                                surf_global_aux_vars(ghosted_id), &
                                option)
    ! [rho*h*T*Cwi]
    xx_loc_p(istart+1) = surf_global_aux_vars(ghosted_id)%den_kg(1)* &
                         xx_loc_p(istart)* &
                         (surf_global_aux_vars(ghosted_id)%temp(1) + 273.15d0)* &
                         surf_th_aux_vars(ghosted_id)%Cwi
  enddo
   
  ! Boundary aux vars
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
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,HET_DIRICHLET)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      surf_global_aux_vars_bc(sum_connection)%temp(1) = xxbc(2)
      call SurfaceTHAuxVarCompute(xxbc, &
                                  surf_th_aux_vars_bc(sum_connection), &
                                  surf_global_aux_vars_bc(sum_connection), &
                                  option)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/Sink aux vars
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

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1

      if (associated(source_sink%flow_condition%temperature)) then
        if(source_sink%flow_condition%temperature%itype/=HET_DIRICHLET) then
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        else
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        endif
      else
        tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+1)
      endif

      xxss = xx_loc_p(istart:iend)
      xxss(2) = tsrc1

      surf_global_aux_vars_ss(sum_connection)%temp(1) = tsrc1
      call SurfaceTHAuxVarCompute(xxss, &
                                  surf_th_aux_vars_ss(sum_connection), &
                                  surf_global_aux_vars_ss(sum_connection), &
                                  option)
    enddo
    source_sink => source_sink%next
  enddo

  patch%surf_aux%SurfaceTH%aux_vars_up_to_date = PETSC_TRUE

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr)

end subroutine SurfaceTHUpdateAuxVars

! ************************************************************************** !
!> This routine updates the temperature after TSSolve.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/25/13
! ************************************************************************** !
subroutine SurfaceTHUpdateTemperature(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module

  implicit none

  type(surface_realization_type) :: surf_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_aux_vars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_aux_vars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: temp
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_global_aux_vars => patch%surf_aux%SurfaceGlobal%aux_vars
  surf_global_aux_vars_bc => patch%surf_aux%SurfaceGlobal%aux_vars_bc
  surf_global_aux_vars_ss => patch%surf_aux%SurfaceGlobal%aux_vars_ss
  surf_aux_vars => patch%surf_aux%SurfaceTH%aux_vars
  surf_aux_vars_bc => patch%surf_aux%SurfaceTH%aux_vars_bc

  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr)

  ! Update internal aux vars
  do ghosted_id = 1,grid%ngmax
    local_id = grid%nG2L(ghosted_id)
    if(local_id>0) then
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      if(xx_loc_p(istart)<1.d-15) then
        temp = 0.d0
      else
        ! T^{t+1} = (rho Cwi hw T)^{t+1} / (rho Cw)^{t} (hw)^{t+1}
        temp = xx_loc_p(iend)/xx_loc_p(istart)/ &
                surf_global_aux_vars(ghosted_id)%den_kg(1)/ &
                surf_aux_vars(ghosted_id)%Cwi - 273.15d0
      endif
      surf_global_aux_vars(ghosted_id)%temp(1) = temp
    endif
  enddo

  ! Update boundary aux vars
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      surf_global_aux_vars_bc(sum_connection)%temp(1) = &
        surf_global_aux_vars(ghosted_id)%temp(1)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Update source/sink aux vars
  source_sink => patch%source_sinks%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      surf_global_aux_vars_ss(sum_connection)%temp(1) = &
        surf_global_aux_vars(ghosted_id)%temp(1)

    enddo
    source_sink => source_sink%next
  enddo

  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr)

end subroutine SurfaceTHUpdateTemperature

! ************************************************************************** !
!> This routine updates solution after a successful time step
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHUpdateSolution(surf_realization)

  use Surface_Realization_class
  use Surface_Field_module

  implicit none

  type(surface_realization_type)   :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode                   :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr)

end subroutine SurfaceTHUpdateSolution


end module Surface_TH_module

#endif
