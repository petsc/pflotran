#ifdef GEOMECH

module Geomechanics_Auxiliary_module

  use Geomechanics_Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: geomech_auxiliary_type
    type(geomech_global_type), pointer :: GeomechGlobal
  end type geomech_auxiliary_type
  
  public :: GeomechAuxInit, &
            GeomechAuxDestroy

contains

! ************************************************************************** !
!
! GeomechAuxInit: Nullifies pointers in geomech auxiliary type
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechAuxInit(geomech_aux)

  implicit none
  
  type(geomech_auxiliary_type) :: geomech_aux
  
  nullify(geomech_aux%GeomechGlobal)
  
end subroutine GeomechAuxInit

! ************************************************************************** !
!
! GeomechAuxDestroy: Strips a geomech auxiliary type
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechAuxDestroy(geomech_aux)

  implicit none
  
  type(geomech_auxiliary_type) :: geomech_aux
  
  call GeomechGlobalAuxDestroy(geomech_aux%GeomechGlobal)

  nullify(geomech_aux%GeomechGlobal)

end subroutine GeomechAuxDestroy

end module Geomechanics_Auxiliary_module

#endif
