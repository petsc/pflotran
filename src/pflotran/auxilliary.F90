module Auxilliary_module
  
  use Richards_Aux_module
  use Richards_Lite_Aux_module
  use Reactive_Transport_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: auxilliary_type 
    type(reactive_transport_type), pointer :: RT
    type(richards_type), pointer :: Richards
    type(richards_lite_type), pointer :: RichardsLite
  end type auxilliary_type
  
  public :: AuxDestroy

contains

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
  call RichardsAuxDestroy(aux%Richards)
  call RichardsLiteAuxDestroy(aux%RichardsLite)
  
  nullify(aux%RT)
  nullify(aux%Richards)
  nullify(aux%RichardsLite)
  
end subroutine AuxDestroy

end module Auxilliary_module
