module Geomechanics_Coupler_module
 
!  use Geomechanics_Condition_module
  use Geomechanics_Region_module
 
  implicit none

  private
 
#include "definitions.h"

  ! coupler types
  ! SK: Note that there is no initial coupler since we solve 
  ! a quasi-static problem for geomechanics (when coupled to flow, otherwise
  ! it is a steady state problem)
  PetscInt, parameter, public :: GM_BOUNDARY_COUPLER_TYPE = 1
  PetscInt, parameter, public :: GM_SRC_SINK_COUPLER_TYPE = 2

   type, public :: geomech_coupler_type
    PetscInt :: id                                         ! id of coupler
    character(len=MAXWORDLENGTH) :: name                   ! name of coupler
    PetscInt :: itype                                      ! integer defining type
    character(len=MAXWORDLENGTH) :: ctype                  ! character string defining type
    character(len=MAXWORDLENGTH) :: geomech_condition_name ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name            ! character string defining name of region to be applied
    PetscInt :: igeomech_condition                         ! id of condition in condition array/list
    PetscInt :: iregion                                    ! id of region in region array/list
    PetscInt, pointer :: geomech_aux_int_var(:,:)          ! auxiliary array for integer value
    PetscReal, pointer :: geomech_aux_real_var(:,:)        ! auxiliary array for real values
 !   type(flow_condition_type), pointer :: flow_condition     ! pointer to condition in condition array/list
    type(gm_region_type), pointer :: region                ! pointer to region in region array/list
    type(geomech_coupler_type), pointer :: next                 ! pointer to next coupler
  end type geomech_coupler_type
  
  type, public :: geomech_coupler_ptr_type
    type(geomech_coupler_type), pointer :: ptr
  end type geomech_coupler_ptr_type
    
  type, public :: geomech_coupler_list_type
    PetscInt :: num_couplers
    type(geomech_coupler_type), pointer :: first
    type(geomech_coupler_type), pointer :: last
    type(geomech_coupler_type), pointer :: array(:)    
  end type geomech_coupler_list_type
  
end module Geomechanics_Coupler_module
