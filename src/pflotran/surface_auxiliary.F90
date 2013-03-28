#ifdef SURFACE_FLOW

module Surface_Auxiliary_module

  use Surface_Global_Aux_module
!  use Surface_Flow_Aux_module
  use Surface_TH_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: surface_auxiliary_type
    type(surface_global_type), pointer :: SurfaceGlobal
    type(surface_th_type), pointer :: SurfaceTH
  end type surface_auxiliary_type
  
  public :: SurfaceAuxInit, &
            SurfaceAuxDestroy

contains

! ************************************************************************** !
!> This routine initializes a surface-auxiliary object
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceAuxInit(surf_aux)

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  nullify(surf_aux%SurfaceGlobal)
  nullify(surf_aux%SurfaceTH)
  
end subroutine SurfaceAuxInit

! ************************************************************************** !
!> This routine deallocates pointers in a surface-auxiliary object
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceAuxDestroy(surf_aux)

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  call SurfaceGlobalAuxDestroy(surf_aux%SurfaceGlobal)
  call SurfaceTHAuxDestroy(surf_aux%SurfaceTH)

  nullify(surf_aux%SurfaceGlobal)
  nullify(surf_aux%SurfaceTH)

end subroutine SurfaceAuxDestroy

end module Surface_Auxiliary_module

#endif
