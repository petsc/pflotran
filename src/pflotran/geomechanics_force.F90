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
! GeomechanicsForceInitialGuess: Sets up the inital guess for the solution
!                                The boudnary conditions are set here
! author: Satish Karra, LANL
! date: 06/19/13
!
! ************************************************************************** !
subroutine GeomechanicsForceInitialGuess(realization)

  use Geomechanics_Realization_module
  use Geomechanics_Field_module
  use Option_module
  
  implicit none
  
  type(geomech_realization_type) :: realization
  
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field
  
  option => realization%option
  field => realization%geomech_field


end subroutine GeomechanicsForceInitialGuess

! ************************************************************************** !
!
! GeomechForceUpdateAuxVars: Updates the geomechanics variables
! author: Satish Karra, LANL
! date: 06/18/13
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
  type(gm_region_type), pointer             :: region
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)

  PetscInt :: ghosted_id, local_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field

  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  
  call GeomechGridVecGetArrayF90(grid,geomech_field%disp_xx_loc,xx_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) = &
        xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
  enddo
   
  call GeomechGridVecRestoreArrayF90(grid,geomech_field%disp_xx_loc, &
                                     xx_loc_p,ierr)

end subroutine GeomechForceUpdateAuxVars

end module Geomechanics_Force_module

#endif
