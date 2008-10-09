module Auxilliary_module
  
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use Mphase_Aux_module
  use Immis_Aux_module
  
  implicit none

  private

#include "definitions.h"

  type, public :: auxilliary_type 
    type(reactive_transport_type), pointer :: RT
    type(thc_type), pointer :: THC
    type(richards_type), pointer :: Richards
    type(mphase_type), pointer :: Mphase
    type(immis_type), pointer :: Immis
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
  nullify(aux%Richards)
  nullify(aux%Mphase)
  nullify(aux%Immis)
  
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
  call RichardsAuxDestroy(aux%Richards)
  !call MphaseAuxDestroy(aux%Mphase)
  
  nullify(aux%RT)
  nullify(aux%THC)
  nullify(aux%Richards)
  nullify(aux%Mphase)
  nullify(aux%Immis)
  
end subroutine AuxDestroy

end module Auxilliary_module
