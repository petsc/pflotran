module Auxiliary_module
  
  use Global_Aux_module
  use TH_Aux_module
  use THC_Aux_module
  use THMC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use Mphase_Aux_module
  use Immis_Aux_module
  use Miscible_Aux_module
  use Flash2_Aux_Module
  use General_Aux_module
  use Material_Aux_module
  use Secondary_Continuum_Aux_module
#ifdef SURFACE_FLOW
  !use Surface_Flow_Aux_module
#endif
  
  implicit none

  private

#include "definitions.h"

  type, public :: auxiliary_type 
    type(global_type), pointer :: Global
    type(reactive_transport_type), pointer :: RT
    type(th_type), pointer :: TH
    type(thc_type), pointer :: THC
    type(thmc_type), pointer :: THMC
    type(richards_type), pointer :: Richards
    type(mphase_type), pointer :: Mphase
    type(immis_type), pointer :: Immis
    type(miscible_type), pointer :: Miscible
    type(flash2_type), pointer :: Flash2
    type(general_type), pointer :: General
    type(material_type), pointer :: Material
    type(sc_heat_type), pointer :: SC_heat
    type(sc_rt_type), pointer :: SC_RT
#ifdef SURFACE_FLOW
    !type(surface_flow_type),pointer :: SurfaceFlow
#endif
  end type auxiliary_type
  
  public :: AuxInit, &
            AuxDestroy

contains

! ************************************************************************** !
!
! AuxInit: Nullifies pointers in auxiliary object
! author: Glenn Hammond
! date: 04/09/08
!
! ************************************************************************** !
subroutine AuxInit(aux)

  implicit none
  
  type(auxiliary_type) :: aux
  
  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%TH)
  nullify(aux%THC)
  nullify(aux%THMC)
  nullify(aux%Richards)
  
  nullify(aux%Mphase)
  nullify(aux%Immis)
  nullify(aux%Flash2)
  nullify(aux%Miscible)
  nullify(aux%General)
  nullify(aux%Material)
  nullify(aux%SC_heat)
  nullify(aux%SC_RT)

#ifdef SURFACE_FLOW
  !nullify(aux%SurfaceFlow)
#endif
end subroutine AuxInit

! ************************************************************************** !
!
! AuxDestroy: Deallocates any allocated pointers in auxiliary object
! author: Glenn Hammond
! date: 04/09/08
!
! ************************************************************************** !
subroutine AuxDestroy(aux)

  implicit none
  
  type(auxiliary_type) :: aux
  
  call GlobalAuxDestroy(aux%Global)
  call RTAuxDestroy(aux%RT)
  call THAuxDestroy(aux%TH)
  call THCAuxDestroy(aux%THC)
  call THMCAuxDestroy(aux%THMC)
  call RichardsAuxDestroy(aux%Richards)
  call MphaseAuxDestroy(aux%Mphase)
  call MiscibleAuxDestroy(aux%Miscible)
  call GeneralAuxDestroy(aux%General)
  call MaterialAuxDestroy(aux%Material)
  call SecondaryAuxHeatDestroy(aux%SC_heat)
  call SecondaryAuxRTDestroy(aux%SC_RT)
  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%THC)
  nullify(aux%THMC)
  nullify(aux%Richards)
  nullify(aux%Mphase)
  nullify(aux%Immis)
  nullify(aux%Miscible)
  nullify(aux%General)
  nullify(aux%Material)
  nullify(aux%SC_Heat)
  nullify(aux%SC_RT)
#ifdef SURFACE_FLOW
  !call SurfaceFlowAuxDestroy(aux%SurfaceFlow)
  !nullify(aux%SurfaceFlow)
#endif
end subroutine AuxDestroy

end module Auxiliary_module
