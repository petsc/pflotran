module Subcontinuum_Transport_module

  use Transport_module
  use Reaction_module
  use Subcontinuum_module
  use Reactive_Transport_Aux_module
  use Subcontinuum_Transport_Aux_module
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
            STSetup, &
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

! ********************************************************************** !
!
! STSetup: 
! author: Jitendra Kumar 
! date: 11/16/2010
!
! Adopted from RTSetup
! ********************************************************************* !
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

! ********************************************************************* !
!
! STSetupPatch: 
! author: Glenn Hammond
! date: 02/22/08
!
! ********************************************************************* !
subroutine STSetupPatch(realization)

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
  PetscInt :: iphase, num_grids, offset
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction
  
  patch%aux%ST%num_cells = patch%nlmax
    
  do iaux=1,patch%aux%ST%num_cells
    patch%aux%ST%st_type1(iaux)%num_subcontinuum =  &
                        patch%num_subcontinuum_type(iaux,1) 
    do jaux=1,patch%aux%ST%st_type1(iaux)%num_subcontinuum
      patch%aux%ST%st_type1(iaux)%st_type(jaux) => RTAuxCreate(option)
      patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%ncomp = reaction%ncomp
      patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%naqcomp = reaction%naqcomp
      patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%nimcomp = 0
      patch%aux%ST%st_type1(iaux)%st_type(jaux)rt_parameter%offset_aq = reaction%offset_aq
      if (reaction%ncollcomp > 0) then
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%ncoll = reaction%ncoll
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%offset_coll = reaction%offset_coll
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%ncollcomp = reaction%ncollcomp
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%offset_collcomp = reaction%offset_collcomp
        allocate(patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%pri_spec_to_coll_spec(reaction%naqcomp))
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%pri_spec_to_coll_spec = &
          reaction%pri_spec_to_coll_spec
        allocate(patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%coll_spec_to_pri_spec(reaction%ncollcomp))
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%coll_spec_to_pri_spec = &
          reaction%coll_spec_to_pri_spec
      endif
    
      ! allocate aux_var data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
      option%iflag = 1 ! allocate mass_balance array
#else  
      option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
      offset = grid%subcontinuum_grid_offset(iaux,2)
      num_subgrids = grid%subcontinuum_grid(offset-1+jaux)
      allocate(patch%aux%ST%st_type1(iaux)%st_type(jaux)%aux_vars(num_subgrids))
      do local_id = 1, num_subgrids 
        call RTAuxVarInit(patch%aux%ST%st_type1(iaux)%st_type(jaux)%aux_vars(local_id),reaction,option)
      enddo
      patch%aux%ST%st_type1(iaux)%st_type(jaux)%num_aux = num_subgrids 
  
      ! count the number of boundary connections and allocate
      ! aux_var data structures for them
      boundary_condition => patch%boundary_conditions%first
      sum_connection = 0    
      !do 
      !  if (.not.associated(boundary_condition)) exit
      !  sum_connection = sum_connection + &
      !               boundary_condition%connection_set%num_connections
      !  boundary_condition => boundary_condition%next
      !enddo
      ! NOTE (Jitu): Hardcoding num of boundary connections to 2 for
      ! subcontinuum problem (Revisit this later.) 11/23/2010
      sum_connection = 2    
  
      if (sum_connection > 0) then
        option%iflag = 1 ! enable allocation of mass_balance array 
        allocate(patch%aux%ST%st_type1(iaux)%st_type(jaux)%aux_vars_bc(sum_connection))
        do iconn = 1, sum_connection
          call RTAuxVarInit(patch%aux%ST%st_type1(iaux)%st_type(jaux)%aux_vars_bc(iconn),reaction,option)
        enddo
      endif
      patch%aux%ST%st_type1(iaux)%st_type(jaux)%num_aux_bc = 2 
      option%iflag = 0

      ! create zero array for zeroing residual and Jacobian (1 on diagonal)
      ! for inactive cells (and isothermal)
      call RTCreateZeroArray(patch,reaction,option)
  
      ! initialize parameters
      cur_fluid_property => realization%fluid_properties
      do 
        if (.not.associated(cur_fluid_property)) exit
        iphase = cur_fluid_property%phase_id
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%diffusion_coefficient(iphase) = &
          cur_fluid_property%diffusion_coefficient
        cur_fluid_property => cur_fluid_property%next
      enddo
  
      if (associated(realization%material_properties)) then
        patch%aux%ST%st_type1(iaux)%st_type(jaux)%rt_parameter%dispersivity = &
          realization%material_properties%longitudinal_dispersivity
      endif
    enddo
  enddo  
end subroutine STSetupPatch


