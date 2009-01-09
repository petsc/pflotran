module Transport_module

  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none
  
  private 

#include "definitions.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"

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
subroutine TFlux(rt_aux_var_up,global_aux_var_up,por_up,tor_up,dist_up, &
                 rt_aux_var_dn,global_aux_var_dn,por_dn,tor_dn,dist_dn, &
                 area,option,velocity,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: Res(option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: sat_up, sat_dn
  PetscReal :: stp_up, stp_dn
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  iphase = 1
  q = velocity(iphase)
  
  sat_up = global_aux_var_up%sat(iphase)
  sat_dn = global_aux_var_dn%sat(iphase)
  
  if (sat_up > eps .and. sat_dn > eps) then
    stp_up = sat_up*tor_up*por_up 
    stp_dn = sat_dn*tor_dn*por_dn
    ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
    weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
    diffusion = weight*(option%disp+option%difaq)
  endif
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)
  
! Add in multiphase, clu 12/29/08
#if 1  
  do
   iphase = iphase +1 
   if (iphase > option%nphase) exit
! super critical CO2 phase have the index 2: need implementation
   q = velocity(iphase)
  
  sat_up = global_aux_var_up%sat(iphase)
  sat_dn = global_aux_var_dn%sat(iphase)
  
  if (sat_up > eps .and. sat_dn > eps) then
    stp_up = sat_up*tor_up*por_up 
    stp_dn = sat_dn*tor_dn*por_dn
    ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
    weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec

   diffusion = 0D0 
   if(iphase ==2) diffusion = weight*(option%disp+option%difsc)

  endif
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  Res(1:option%ntrandof) = Res (1:option%ntrandof) + & 
                      coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                      coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)
 enddo
#endif

end subroutine TFlux

! ************************************************************************** !
!
! TFluxDerivative: Computes derivatives of flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFluxDerivative(rt_aux_var_up,global_aux_var_up,por_up,tor_up,dist_up, &
                           rt_aux_var_dn,global_aux_var_dn,por_dn,tor_dn,dist_dn, &
                           area,option,velocity,J_up,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: J_up(option%ntrandof,option%ntrandof), J_dn(option%ntrandof,option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: sat_up, sat_dn
  PetscReal :: stp_up, stp_dn
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  
  diffusion = 0.d0

  iphase = 1
  q = velocity(iphase)

  sat_up = global_aux_var_up%sat(iphase)
  sat_dn = global_aux_var_dn%sat(iphase)
    
  if (sat_up > eps .and. sat_dn > eps) then
    stp_up = sat_up*tor_up*por_up
    stp_dn = sat_dn*tor_dn*por_dn
    ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
    weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
    diffusion = weight*(option%disp+option%difaq)
  endif
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)
  !       = m^3 water/sec
  coef_up = coef_up*area
  coef_dn = coef_dn*area

  ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  if (associated(rt_aux_var_dn%dtotal)) then
    J_up = rt_aux_var_up%dtotal(:,:,iphase)*coef_up*1000.d0
    J_dn = rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
  else  
    J_up = 0.d0
    J_dn = 0.d0
    do icomp = 1, option%ntrandof
      J_up(icomp,icomp) = coef_up*global_aux_var_up%den_kg(iphase)
      J_dn(icomp,icomp) = coef_dn*global_aux_var_dn%den_kg(iphase)
    enddo
  endif

! Add in multiphase, clu 12/29/08
#if 1  
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit
! super critical CO2 phase
    q = velocity(iphase)

    sat_up = global_aux_var_up%sat(iphase)
    sat_dn = global_aux_var_dn%sat(iphase)
    
    if (sat_up > eps .and. sat_dn > eps) then
      stp_up = sat_up*tor_up*por_up
      stp_dn = sat_dn*tor_dn*por_dn
      ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
      weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
      if(iphase==2) diffusion = weight*(option%disp+option%difsc)
    endif
  
    !upstream weighting
    ! units = (m^3 water/m^2 bulk/sec)
    if (q > 0.d0) then
      coef_up =  diffusion+q
      coef_dn = -diffusion
    else
      coef_up =  diffusion
      coef_dn = -diffusion+q
    endif
  
    ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)
    !       = m^3 water/sec
    coef_up = coef_up*area
    coef_dn = coef_dn*area

    ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
    if (associated(rt_aux_var_dn%dtotal)) then
      J_up = J_up + rt_aux_var_up%dtotal(:,:,iphase)*coef_up*1000.d0
      J_dn = J_dn + rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
    else  
   !   J_up = 0.d0
   !   J_dn = 0.d0
      do icomp = 1, option%ntrandof
        J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up*global_aux_var_up%den_kg(iphase)
        J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
      enddo
    endif
  enddo
#endif

end subroutine TFluxDerivative

! ************************************************************************** !
!
! TBCFlux: Computes boundary flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TBCFlux(ibndtype, &
                   rt_aux_var_up,global_aux_var_up, &
                   rt_aux_var_dn,global_aux_var_dn, &
                   por_dn,tor_dn,dist_dn, &
                   area,option,velocity,Res)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_dn, tor_dn,  dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: Res(option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: coef_up, coef_dn
  PetscReal :: diffusion, q
  PetscReal :: sat_up, sat_dn
  
  diffusion = 0.d0

  iphase = 1
  q = velocity(iphase)
  
  sat_up = global_aux_var_up%sat(iphase)
  sat_dn = global_aux_var_dn%sat(iphase)

  select case(ibndtype)
    case(DIRICHLET_BC)
      if (sat_up > eps .and. sat_dn > eps) then
        ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
        weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
        ! need to account for multiple phases
        ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
        diffusion = weight*(option%disp+option%difaq)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        ! same as dirichlet above
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
          diffusion = weight*(option%disp+option%difaq)
        endif    
      endif
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif

  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0

  ! units = (L water/sec)*(mol/L) = mol/s  
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)  

! Add in multiphase, clu 12/29/08

#if 1  
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit
      q = velocity(iphase)
  
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)

      select case(ibndtype)
        case(DIRICHLET_BC)
          if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
            weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
           ! need to account for multiple phases
           ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
           if( iphase == 2) diffusion = weight*(option%disp+option%difaq)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
           ! same as dirichlet above
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
            diffusion =0D0
            if( iphase == 2) diffusion = weight*(option%disp+option%difaq)
            endif    
          endif
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion+q
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion+q
  endif

  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0

  ! units = (L water/sec)*(mol/L) = mol/s  
  Res(1:option%ntrandof) = Res(1:option%ntrandof) + &
                           coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)  
 enddo
#endif
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
                             rt_aux_var_up,global_aux_var_up, &
                             rt_aux_var_dn,global_aux_var_dn, &
                             por_dn,tor_dn,dist_dn, &
                             area,option,velocity,J_dn)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  PetscReal :: J_dn(option%ntrandof,option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: coef_dn
  PetscReal :: diffusion, q
  PetscReal :: sat_up, sat_dn  
  
  diffusion = 0.d0

  iphase = 1
  q = velocity(iphase)
  
  sat_up = global_aux_var_up%sat(iphase)
  sat_dn = global_aux_var_dn%sat(iphase)
  
  select case(ibndtype)
    case(DIRICHLET_BC)
      if (sat_up > eps .and. sat_dn > eps) then
        ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
        weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
        ! need to account for multiple phases
        ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
        diffusion = weight*(option%disp+option%difaq)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
          diffusion = weight*(option%disp+option%difaq)
        endif  
      endif  
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)  
  if (q > 0.d0) then
    coef_dn = -diffusion
  else
    coef_dn = -diffusion+q
  endif

  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)
  !       = m^3 water/sec  
  coef_dn = coef_dn*area

  ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  if (associated(rt_aux_var_dn%dtotal)) then
    J_dn = rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
  else
    J_dn = 0.d0
    do icomp = 1, option%ntrandof
      J_dn(icomp,icomp) = coef_dn*global_aux_var_dn%den_kg(iphase)
    enddo
  endif

! Add in multiphase, clu 12/29/08
#if 1  
  do 
    iphase = iphase + 1
    if (iphase > option%nphase) exit

! super critical CO2 phase
      q = velocity(iphase)
  
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)
  
      select case(ibndtype)
        case(DIRICHLET_BC)
          if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
            weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
           ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
            diffusion = weight*(option%disp+option%difaq)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
            ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
              diffusion = weight*(option%disp+option%difaq)
           endif  
          endif  
        case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
      end select

      !upstream weighting
      ! units = (m^3 water/m^2 bulk/sec)  
      if (q > 0.d0) then
        coef_dn = -diffusion
      else
        coef_dn = -diffusion+q
      endif

      ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)
      !       = m^3 water/sec  
      coef_dn = coef_dn*area

      ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
      if (associated(rt_aux_var_dn%dtotal)) then
        J_dn = rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
      else
        J_dn = 0.d0
        do icomp = 1, option%ntrandof
          J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
        enddo
       endif
   enddo
#endif

end subroutine TBCFluxDerivative

end module Transport_module
