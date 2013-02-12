#if 0
!#ifdef SURFACE_FLOW

module Surface_Flow_Aux_module

  implicit none
  
  private
  
#include "definitions.h"
  PetscInt, parameter, public :: KINEMATIC_WAVE = 1
  PetscInt, parameter, public :: DIFFUSION_WAVE = 2
  
  type, public :: surface_flow_type
    PetscInt :: formulation_type
  end type surface_flow_type
  
  public :: SurfaceFlowAuxCreate, &
            SurfaceFlowAuxDestroy
  
contains

! ************************************************************************** !
!> This function allocates and initializes a surfaceflow auxilliary object
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
function SurfaceFlowAuxCreate()

  implicit none
  
  type(surface_flow_type), pointer :: SurfaceFlowAuxCreate
  type(surface_flow_type), pointer :: aux
  
  aux%formulation_type = KINEMATIC_WAVE
  
  SurfaceFlowAuxCreate => aux
  
end function SurfaceFlowAuxCreate

! ************************************************************************** !
!> This routine deallocates a surfaceflow auxilliary object
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 08/03/12
! ************************************************************************** !
subroutine SurfaceFlowAuxDestroy(aux)

  implicit none
  
  type(surface_flow_type), pointer :: aux

  if(.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)

end subroutine SurfaceFlowAuxDestroy

end module Surface_Flow_Aux_module

#endif
