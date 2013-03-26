#ifdef SURFACE_FLOW

module Surface_TH_module

  use Surface_Global_Aux_module
  use Surface_TH_Aux_module
  
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
         SurfaceTHUpdateSolution

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
  use Level_module
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

  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,NFLOWDOF)
  
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

        call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)

        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          iend = ghosted_id*option%nflowdof
          istart = iend-option%nflowdof+1
          vec_p((iconn-1)*option%nflowdof+1:iconn*option%nflowdof) = &
            xx_loc_p(istart:iend)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)
        call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)

        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_ndof, &
                        surf_field%subsurf_temp_vec_ndof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_ndof, &
                        surf_field%subsurf_temp_vec_ndof,surf_field%press_subsurf, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)

      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    endif
    coupler => coupler%next
  enddo

end subroutine SurfaceTHUpdateSurfBC

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
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class

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
                             surf_field%vol_subsurf_2_surf, &
                             surf_field%subsurf_temp_vec_ndof, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_ndof, &
                           surf_field%vol_subsurf_2_surf, &
                           surf_field%subsurf_temp_vec_ndof, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)

        call VecGetArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)
        do iconn=1,coupler%connection_set%num_connections
          ! Flux of water
          coupler%flow_aux_real_var(ONE_INTEGER,iconn)=-vec_p((iconn-1)*option%nflowdof+1)/dt*den
          
          ! Heat flux
          coupler%flow_aux_real_var(TWO_INTEGER,iconn)=-vec_p(iconn*option%nflowdof)/dt*den
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_ndof,vec_p,ierr)

        call VecSet(surf_field%vol_subsurf_2_surf,0.d0,ierr)
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
  use Level_module
  use Patch_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Surface_Field_module
  use Water_EOS_module
  use Discretization_module
  use Surface_Realization_class
  use Realization_Base_class

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
  use Level_module
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
  
  Vec            :: destin_mpi_vec, source_mpi_vec
  PetscErrorCode :: ierr
  PetscReal :: den          ! density      [kg/m^3]
  PetscInt :: local_id,iconn
  
  PetscReal, pointer :: hw_p(:)   ! head [m]
  PetscReal, pointer :: press_sub_p(:) ! Pressure [Pa]
  PetscReal, pointer :: sat_func_id_p(:)
  PetscReal, pointer :: Dq_p(:)
  PetscReal, pointer :: vol_p(:)
  PetscReal, pointer :: area_p(:)
  PetscReal, pointer :: dist_p(:)
  
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
    
  PetscBool :: coupler_found = PETSC_FALSE
  PetscBool :: v_darcy_limit

  patch      => realization%patch
  surf_patch => surf_realization%patch
  option     => realization%option
  grid       => realization%discretization%grid
  field      => realization%field
  surf_grid  => surf_realization%discretization%grid
  surf_field => surf_realization%surf_field

  call density(option%reference_temperature,option%reference_pressure,den)

  call GridVecGetArrayF90(surf_grid,surf_field%press_subsurf,press_sub_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%flow_xx_loc,hw_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%sat_func_id,sat_func_id_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%Dq,Dq_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%vol_subsurf_2_surf,vol_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

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
      
      do local_id=1,surf_grid%nlmax
        press_surf=hw_p(local_id)*(abs(option%gravity(3)))*den+option%reference_pressure

        press_up = press_sub_p(local_id)
        press_dn = press_surf
        gravity = dist_p(local_id)*den
        
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
          patch%saturation_function_array(int(sat_func_id_p(local_id)))%ptr, &
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
          ! Flow is happening from surface to subsurface
          if ( abs(v_darcy) > hw_p(local_id)/option%surf_flow_dt ) then
            v_darcy = -hw_p(local_id)/option%surf_flow_dt
            v_darcy_limit=PETSC_TRUE
          endif
        else
          ! Exfiltration is occuring
        endif
        
        vol_p(local_id)=vol_p(local_id)+v_darcy*area_p(local_id)*option%surf_flow_dt
        coupler%flow_aux_real_var(ONE_INTEGER,local_id)=v_darcy
        if(abs(v_darcy)>v_darcy_max) v_darcy_max=v_darcy
      enddo

    endif
    
    coupler => coupler%next
  enddo
  
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%vol_subsurf_2_surf,vol_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%Dq,Dq_p,ierr)  
  call GridVecRestoreArrayF90(surf_grid,surf_field%sat_func_id,sat_func_id_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%flow_xx_loc,hw_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%press_subsurf,press_sub_p,ierr)

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
  use Level_module
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
  PetscReal, pointer :: dist_p(:)
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
        call GridVecGetArrayF90(grid,field%perm_xx_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%perm_xx_loc,xx_loc_p, ierr)
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
        call GridVecGetArrayF90(grid,field%perm_yy_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%perm_yy_loc,xx_loc_p, ierr)
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
        call GridVecGetArrayF90(grid,field%perm_zz_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%perm_zz_loc,xx_loc_p, ierr)
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
        call GridVecGetArrayF90(grid,field%porosity_loc,xx_loc_p, ierr)
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=xx_loc_p(ghosted_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%porosity_loc,xx_loc_p, ierr)
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

        ! sat_func_id
        call VecGetArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          vec_p(iconn)=patch%sat_func_id(ghosted_id)
        enddo
        call VecRestoreArrayF90(surf_field%subsurf_temp_vec_1dof,vec_p,ierr)
        ! Scatter the data
        call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                            surf_field%subsurf_temp_vec_1dof, &
                            surf_field%sat_func_id, &
                            INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                           surf_field%subsurf_temp_vec_1dof, &
                           surf_field%sat_func_id, &
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

  call GridVecGetArrayF90(surf_grid,surf_field%subsurf_xx,xx_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%subsurf_yy,yy_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%subsurf_zz,zz_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%perm_xx,perm_xx_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%perm_yy,perm_yy_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%perm_zz,perm_zz_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%Dq,Dq_p,ierr)
  call GridVecGetArrayF90(surf_grid,surf_field%surf2subsurf_dist_gravity,dist_p,ierr)

  do local_id=1,surf_grid%nlmax
    dist_x = (xx_p(local_id) - surf_grid%x(local_id))
    dist_y = (yy_p(local_id) - surf_grid%y(local_id))
    dist_z = (zz_p(local_id) - surf_grid%z(local_id))
      
    dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z)

    Dq_p(local_id) = (perm_xx_p(local_id)*abs(dist_x)/dist + &
                      perm_yy_p(local_id)*abs(dist_y)/dist + &
                      perm_zz_p(local_id)*abs(dist_z)/dist)/dist
    dist_p(local_id) = dist*(dist_x*option%gravity(1)+ &
                             dist_y*option%gravity(2)+ &
                             dist_z*option%gravity(3))
  enddo

  call GridVecRestoreArrayF90(surf_grid,surf_field%surf2subsurf_dist_gravity,dist_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%subsurf_xx,xx_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%subsurf_yy,yy_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%subsurf_zz,zz_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%perm_xx,perm_xx_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%perm_yy,perm_yy_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%perm_zz,perm_zz_p,ierr)
  call GridVecRestoreArrayF90(surf_grid,surf_field%Dq,Dq_p,ierr)

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

  use water_eos_module
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
  PetscReal :: max_allowable_dt
  PetscReal :: qsrc, qsrc_flow
  PetscReal :: tsrc

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
                                   surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)
  ! Then, update the aux vars
  call SurfaceTHUpdateAuxVars(surf_realization)
  ! override flags since they will soon be out of date  
  patch%surf_aux%SurfaceTH%aux_vars_up_to_date = PETSC_FALSE

  call GridVecGetArrayF90(grid,ff,ff_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%mannings_loc,mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

  ff_p = 0.d0
  Res  = 0.d0
  max_allowable_dt = 1.d10

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
#ifdef STORE_FLOWRATE
      patch%surf_internal_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(TH_PRESSURE_DOF)
      !patch%surf_internal_fluxes(TH_TEMPERATURE_DOF,sum_connection) =
#endif
      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)

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

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         surf_aux_vars_bc(sum_connection), &
                         surf_global_aux_vars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
#ifdef STORE_FLOWRATE
      patch%surf_boundary_fluxes(TH_PRESSURE_DOF,sum_connection) = Res(1)
#endif
      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)
      
      iend = local_id_dn*option%nflowdof
      istart = iend-option%nflowdof+1
      ff_p(istart:iend) = ff_p(istart:iend) + Res(1)/area_p(local_id_dn)
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
    qsrc_flow = source_sink%flow_condition%rate%flow_dataset%time_series%cur_value(1)
      
    if(source_sink%flow_condition%temperature%itype/=HET_DIRICHLET) &
      tsrc = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)

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
      
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ff_p(istart) = ff_p(istart) + qsrc/area_p(local_id)
      ff_p(iend) = 0.d0

    enddo
    source_sink => source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,ff,ff_p, ierr)
  call GridVecRestoreArrayF90(grid,surf_field%mannings_loc,mannings_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)


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

  use water_eos_module
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

  call GridVecGetArrayF90(grid,surf_field%mannings_loc,mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid,surf_field%area,area_p,ierr)

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

      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)

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

      if(abs(vel)>eps) max_allowable_dt = min(max_allowable_dt,dist/abs(vel)/2.d0)
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  call GridVecRestoreArrayF90(grid,surf_field%mannings_loc,mannings_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,surf_field%area,area_p,ierr)
  
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
  PetscReal :: dhead

  ! initialize
  flux = 0.d0

  head_up = surf_global_aux_var_up%head(1) + zc_up
  head_dn = surf_global_aux_var_dn%head(1) + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (surf_global_aux_var_up%head(1)>eps) then
      hw_half = surf_global_aux_var_up%head(1)
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
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
    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)
  endif

  flux = hw_half*vel
  Res(TH_PRESSURE_DOF) = flux*length
  
  Res(TH_TEMPERATURE_DOF) = 0.d0

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
  
  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_aux_var%head(1)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
      else
        vel = -sqrt(dabs(slope))/mannings*((head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = head*vel
  Res(TH_PRESSURE_DOF) = flux*length

  Res(TH_TEMPERATURE_DOF) = 0.d0

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
  
  call GridVecGetArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p, ierr)

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

      if(source_sink%flow_condition%temperature%itype/=HET_DIRICHLET) then
        tsrc1 = source_sink%flow_condition%temperature%flow_dataset%time_series%cur_value(1)
      else
        tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      endif

      xxss = xx_loc_p(istart:iend)
      xxss(2) = tsrc1
      call SurfaceTHAuxVarCompute(xxss, &
                                  surf_th_aux_vars_ss(sum_connection), &
                                  surf_global_aux_vars_ss(sum_connection), &
                                  option)

    enddo
    source_sink => source_sink%next
  enddo

  patch%surf_aux%SurfaceTH%aux_vars_up_to_date = PETSC_TRUE

  call GridVecRestoreArrayF90(grid,surf_field%flow_xx_loc,xx_loc_p, ierr)

end subroutine SurfaceTHUpdateAuxVars

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