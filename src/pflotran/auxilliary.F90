module Auxilliary_module
  
  use THC_Aux_module
  use Richards_Lite_Aux_module
  use Reactive_Transport_Aux_module
  use Mphase_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: auxilliary_type 
    type(reactive_transport_type), pointer :: RT
    type(thc_type), pointer :: THC
    type(richards_lite_type), pointer :: RichardsLite
    type(mphase_type), pointer :: Mphase
  end type auxilliary_type
  
  public :: AuxInit, &
            AuxDestroy

contains

! ************************************************************************** !
!
! AuxInit: Nullifies pointers in auxilliary object
! author: Glenn Hammond
! date: 04/09/08
!
! ************************************************************************** !
subroutine AuxInit(aux)

  implicit none
  
  type(auxilliary_type) :: aux
  
  nullify(aux%RT)
  nullify(aux%THC)
  nullify(aux%RichardsLite)
  nullify(aux%Mphase)
  
end subroutine AuxInit

! ************************************************************************** !
!
! AuxDestroy: Deallocates any allocated pointers in auxilliary object
! author: Glenn Hammond
! date: 04/09/08
!
! ************************************************************************** !
subroutine AuxDestroy(aux)

  implicit none
  
  type(auxilliary_type) :: aux
  
  call RTAuxDestroy(aux%RT)
  call THCAuxDestroy(aux%THC)
  call RichardsLiteAuxDestroy(aux%RichardsLite)
  !call MphaseAuxDestroy(aux%Mphase)
  
  nullify(aux%RT)
  nullify(aux%THC)
  nullify(aux%RichardsLite)
  nullify(aux%Mphase)
  
end subroutine AuxDestroy

end module Auxilliary_module
