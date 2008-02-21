module Transport_module

  implicit none
  
  private 

#include "definitions.h"
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petsclog.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  
  public :: TFlux, TFluxDerivative, TBCFlux, TBCFluxDerivative
  
contains

! ************************************************************************** !
!
! TFlux: Computes flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFlux(aux_var_up,por_up,tor_up,sat_up,dist_up, &
                 aux_var_dn,por_dn,tor_dn,sat_dn,dist_dn, &
                 area,option,velocity,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var_up, aux_var_dn
  PetscReal :: por_up, tor_up, sat_up, dist_up
  PetscReal :: por_dn, tor_dn, sat_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: weight
  PetscReal :: stp_up, stp_dn
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  q = velocity(1)
  if (sat_up > eps .and. sat_dn > eps) then
    stp_up = sat_up*tor_up*por_up
    stp_dn = sat_dn*tor_dn*por_dn
    weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    diffusion = weight*(option%disp+option%difaq)
  endif
  
  !upstream weighting
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif
  
  coef_up = coef_up*area
  coef_dn = coef_dn*area
  
  Res(1:option%ncomp) = coef_up*aux_var_up%total(1:option%ncomp) + &
                        coef_dn*aux_var_dn%total(1:option%ncomp)
                        
end subroutine TFlux

! ************************************************************************** !
!
! TFluxDerivative: Computes derivatives of flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFluxDerivative(aux_var_up,por_up,tor_up,sat_up,dist_up, &
                           aux_var_dn,por_dn,tor_dn,sat_dn,dist_dn, &
                           area,option,velocity,J_up,J_dn)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: aux_var_up, aux_var_dn
  PetscReal :: por_up, tor_up, sat_up, dist_up
  PetscReal :: por_dn, tor_dn, sat_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: J_up(option%ncomp,option%ncomp), J_dn(option%ncomp,option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: weight
  PetscReal :: stp_up, stp_dn
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  q = velocity(1)
  if (sat_up > eps .and. sat_dn > eps) then
    stp_up = sat_up*tor_up*por_up
    stp_dn = sat_dn*tor_dn*por_dn
    weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    diffusion = weight*(option%disp+option%difaq)
  endif
  
  !upstream weighting
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif
  
  coef_up = coef_up*area
  coef_dn = coef_dn*area
  
  do icomp = 1, option%ncomp
    J_up(icomp,icomp) = coef_up
    J_dn(icomp,icomp) = coef_dn
  enddo

end subroutine TFluxDerivative

! ************************************************************************** !
!
! TBCFlux: Computes boundary flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TBCFlux(ibndtype, &
                   aux_var_up,aux_var_dn,por_dn,tor_dn,sat_dn,dist_dn, &
                   area,option,velocity,Res)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: aux_var_up, aux_var_dn
  PetscReal :: por_dn, tor_dn, sat_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: Res(option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: weight
  PetscReal :: stp_dn
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  q = velocity(1)

  select case(ibndtype)
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      if (sat_dn > eps) then
        stp_dn = sat_dn*tor_dn*por_dn
        weight = stp_dn/dist_dn
        ! need to account for multiple phases
        diffusion = weight*(option%disp+option%difaq)
      endif
    case(NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif

  coef_up = coef_up*area
  coef_dn = coef_dn*area
  
  Res(1:option%ncomp) = coef_up*aux_var_up%total(1:option%ncomp) + &
                        coef_dn*aux_var_dn%total(1:option%ncomp)  

end subroutine TBCFlux

! ************************************************************************** !
!
! TBCFluxDerivative: Computes derivative of boundary flux term in residual 
!                    function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TBCFluxDerivative(ibndtype, &
                             aux_var_up, &
                             aux_var_dn,por_dn,tor_dn,sat_dn,dist_dn, &
                             area,option,velocity,J_dn)

  use Reactive_Transport_Aux_module
  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: aux_var_up, aux_var_dn
  PetscReal :: por_dn, tor_dn, sat_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: J_dn(option%ncomp,option%ncomp)
  
  PetscInt :: icomp
  PetscReal :: weight
  PetscReal :: stp_dn
  PetscReal :: coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  q = velocity(1)

  select case(ibndtype)
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      if (sat_dn > eps) then
        stp_dn = sat_dn*tor_dn*por_dn
        weight = stp_dn/dist_dn
        ! need to account for multiple phases
        diffusion = weight*(option%disp+option%difaq)
      endif
    case(NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  if (q > 0.d0) then
    coef_dn = -diffusion
  else
    coef_dn = -diffusion+q
  endif
  
  coef_dn = coef_dn*area

  do icomp = 1, option%ncomp
    J_dn(icomp,icomp) = coef_dn
  enddo

end subroutine TBCFluxDerivative

end module Transport_module
