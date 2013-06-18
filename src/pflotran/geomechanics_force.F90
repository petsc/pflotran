#ifdef GEOMECH 

module Geomechanics_Force_module

  use Geomechanics_Global_Aux_module
  
  implicit none
  
  private
  
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
!#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"
#include "finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: GeomechForceSetup, &
            GeomechForceUpdateAuxVars
  
contains

! ************************************************************************** !
!
! GeomechForceSetup: Sets up the geomechanics calculations
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechForceSetup(geomech_realization)

  use Geomechanics_Realization_module
  
  type(geomech_realization_type) :: geomech_realization

  call GeomechForceSetPlotVariables(geomech_realization)
  
end subroutine GeomechForceSetup

! ************************************************************************** !
!
! GeomechForceSetPlotVariables: Set up of geomechanics plot variables
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechForceSetPlotVariables(geomech_realization)
  
  use Geomechanics_Realization_module
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => geomech_realization%output_option%output_variable_list
  
  name = 'disp_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_X)
                               
  name = 'disp_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_Y)
                               
  name = 'disp_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GEOMECH_DISP_Z)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine GeomechForceSetPlotVariables

! ************************************************************************** !
!
! GeomechForceSetPlotVariables: Updates the geomechanics variables
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechForceUpdateAuxVars(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Field_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Material_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Region_module

  implicit none

  type(geomech_realization_type)            :: geomech_realization
  
  type(option_type), pointer                :: option
  type(geomech_patch_type), pointer         :: patch
  type(geomech_grid_type), pointer          :: grid
  type(geomech_field_type), pointer         :: geomech_field
  type(geomech_coupler_type), pointer       :: boundary_condition
  type(geomech_coupler_type), pointer       :: source_sink
  type(gm_region_type), pointer             :: region
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars_bc(:)
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars_ss(:)

  PetscInt :: ghosted_id, local_id, total_verts, ivertex
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(3)
  PetscReal :: xxss(3)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field

  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  geomech_global_aux_vars_bc => patch%geomech_aux%GeomechGlobal%aux_vars_bc
  geomech_global_aux_vars_ss => patch%geomech_aux%GeomechGlobal%aux_vars_ss
  
  call GeomechGridVecGetArrayF90(grid,geomech_field%disp_xx_loc,xx_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(ONE_INTEGER) = &
      xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(TWO_INTEGER) = &
      xx_loc_p(TWO_INTEGER + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(THREE_INTEGER) = &
      xx_loc_p(THREE_INTEGER + (ghosted_id-1)*THREE_INTEGER)
  enddo
   
  ! Boundary aux vars
  boundary_condition => patch%geomech_boundary_conditions%first
  total_verts = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            geomech_global_aux_vars_bc(total_verts)%disp_x = boundary_condition% &
              geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(1) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER) ! sk: check this
        end select
      endif
      
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            xxbc(2) = boundary_condition% &
              geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(2) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER) ! sk: check this
        end select
      endif

      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            xxbc(3) = boundary_condition% &
              geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(3) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER) ! sk: check this
        end select
      endif
      
      geomech_global_aux_vars_bc(total_verts)%disp_vector(:) = xxbc(:)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/Sink aux vars
  ! source/sinks
  source_sink => patch%geomech_source_sinks%first
  total_verts = 0
  do
    if (.not.associated(source_sink)) exit
    region => source_sink%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      xxss(ONE_INTEGER) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER)
      xxss(TWO_INTEGER) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*TWO_INTEGER)
      xxss(THREE_INTEGER) = xx_loc_p(ONE_INTEGER + (ghosted_id-1)*THREE_INTEGER)

      geomech_global_aux_vars_ss(total_verts)%disp_vector(:) = xxss(:)
    enddo
    source_sink => source_sink%next
  enddo

  call GeomechGridVecRestoreArrayF90(grid,geomech_field%disp_xx_loc, &
                                     xx_loc_p,ierr)

end subroutine GeomechForceUpdateAuxVars

end module Geomechanics_Force_module

#endif
