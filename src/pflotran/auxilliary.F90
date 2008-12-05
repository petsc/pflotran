module Auxilliary_module
  
  use Global_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use Mphase_Aux_module
  use Immis_Aux_module
  
  implicit none

  private

#include "definitions.h"

  type, public :: auxilliary_type 
    type(global_type), pointer :: Global
    type(reactive_transport_type), pointer :: RT
    type(thc_type), pointer :: THC
    type(richards_type), pointer :: Richards
    type(mphase_type), pointer :: Mphase
    type(immis_type), pointer :: Immis
  end type auxilliary_type
  
#if 0
  type, public :: auxilliary_coupler_type 
    type(global_auxvar_type), pointer :: global_auxvar
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  end type auxilliary_coupler_type
#endif  
  
  public :: AuxInit, &
            AuxDestroy
!            AuxCouplerInit, &
!            AuxCouplerDestroy

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
  
  nullify(aux%Global)
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
  
  call GlobalAuxDestroy(aux%Global)
  call RTAuxDestroy(aux%RT)
  call THCAuxDestroy(aux%THC)
  call RichardsAuxDestroy(aux%Richards)
  !call MphaseAuxDestroy(aux%Mphase)
  
  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%THC)
  nullify(aux%Richards)
  nullify(aux%Mphase)
  nullify(aux%Immis)
  
end subroutine AuxDestroy

#if 0
! ************************************************************************** !
!
! AuxCouplerInit: Nullifies pointers in auxilliary object for a coupler
! author: Glenn Hammond
! date: 04/09/08
!
! ************************************************************************** !
subroutine AuxCouplerInit(aux)

  implicit none
  
  type(auxilliary_coupler_type) :: aux
  
  nullify(aux%global_auxvar)
  nullify(aux%rt_auxvar)
  
end subroutine AuxCouplerInit

! ************************************************************************** !
!
! AuxCouplerDestroy: Deallocates any allocated pointers in auxilliary object
!                    for a coupler
! author: Glenn Hammond
! date: 12/05/08
!
! ************************************************************************** !
subroutine AuxCouplerDestroy(aux)

  implicit none
  
  type(auxilliary_coupler_type) :: aux
  
  call GlobalAuxVarDestroy(aux%global_auxvar)
  call RTAuxVarDestroy(aux%rt_auxvar)
  
  nullify(aux%global_auxvar)
  nullify(aux%rt_auxvar)
  
end subroutine AuxCouplerDestroy
#endif

end module Auxilliary_module
