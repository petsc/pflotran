module Transport_module

  use Reactive_Transport_Aux_module
  use Global_Aux_module
#ifdef REVISED_TRANSPORT
  use Matrix_Block_Aux_module  
#endif  

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
  
  public :: TDiffusion, &
            TDiffusionBC, &
            TFlux, &
            TFluxDerivative, &
#ifndef REVISED_TRANSPORT
            TBCFlux, &
            TBCFluxDerivative, &
            TFluxAdv, &
            TFluxDerivativeAdv, &
            TBCFluxAdv, &
            TBCFluxDerivativeAdv, &
            TFluxDiff, &
            TFluxDerivativeDiff, &
            TBCFluxDiff, &
            TBCFluxDerivativeDiff, &
#endif
            TFluxCoef
  
contains

! ************************************************************************** !
!
! TDiffusion: Computes diffusion term at cell interface
! author: Glenn Hammond
! date: 02/24/10
!
! ************************************************************************** !
subroutine TDiffusion(global_aux_var_up,por_up,tor_up,dist_up, &
                      global_aux_var_dn,por_dn,tor_dn,dist_dn, &
                      rt_parameter,option,velocity,diffusion)

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: velocity(*)
  type(option_type) :: option
  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: diffusion(option%nphase)
  
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: sat_up, sat_dn
  PetscReal :: stp_up, stp_dn
  PetscReal :: q
  
  diffusion(:) = 0.d0

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
<<<<<<< local
    diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
    diffusion(iphase) = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                        weight*rt_parameter%diffusion_coefficient(iphase)
  endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
        if(iphase ==2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase ==2) diffusion(iphase) = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                               weight*rt_parameter%diffusion_coefficient(iphase)
      endif
    enddo
  endif
#endif  
  
end subroutine TDiffusion

! ************************************************************************** !
!
! TDiffusionBC: Computes diffusion term at cell boundary interface
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TDiffusionBC(ibndtype,global_aux_var_up,global_aux_var_dn, &
                        por_dn,tor_dn,dist_dn, &
                        rt_parameter,option,velocity,diffusion)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_dn, tor_dn,  dist_dn
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  PetscReal :: diffusion(option%nphase)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: weight
  PetscReal :: q
  PetscReal :: sat_up, sat_dn
  
  diffusion(:) = 0.d0

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
<<<<<<< local
        diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
        diffusion(iphase) = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                            weight*rt_parameter%diffusion_coefficient(iphase)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        ! same as dirichlet above
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
          diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
          diffusion(iphase) = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                              weight*rt_parameter%diffusion_coefficient(iphase)
        endif    
      endif
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
            if( iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            if( iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                       weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
         ! same as dirichlet above
            if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
        ! need to account for multiple phases
        ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
              if(iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              if(iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                        weight*rt_parameter%diffusion_coefficient(iphase)
            endif    
          endif
        case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
      end select
    enddo
  endif
#endif

end subroutine TDiffusionBC
#ifndef REVISED_TRANSPORT
! ************************************************************************** !
!
! TFlux: Computes flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFlux(rt_aux_var_up,global_aux_var_up,por_up,tor_up,dist_up, &
                 rt_aux_var_dn,global_aux_var_dn,por_dn,tor_dn,dist_dn, &
                 area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
    diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
    diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                weight*rt_parameter%diffusion_coefficient(iphase)
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do
      iphase = iphase +1 
      if (iphase > option%nphase) exit
! super critical CO2 phase have the index 2: need implementation
      q = velocity(iphase)
      diffusion = 0.d0
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)
  
      if (sat_up > eps .and. sat_dn > eps) then
        stp_up = sat_up*tor_up*por_up 
        stp_dn = sat_dn*tor_dn*por_dn
    ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
        weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
        if(iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                              weight*rt_parameter%diffusion_coefficient(iphase)
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
  endif
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
                           area,rt_parameter,option,velocity,J_up,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
    diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
    diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                weight*rt_parameter%diffusion_coefficient(iphase)

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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit
! super critical CO2 phase
      q = velocity(iphase)
      diffusion = 0D0
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)
    
      if (sat_up > eps .and. sat_dn > eps) then
        stp_up = sat_up*tor_up*por_up
        stp_dn = sat_dn*tor_dn*por_dn
      ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
        weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
        if(iphase==2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase==2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                                weight*rt_parameter%diffusion_coefficient(iphase)
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
        print *,'Dtotal needed for SC problem. STOP'
        stop 
   !   J_up = 0.d0
   !   J_dn = 0.d0
   !   do icomp = 1, option%ntrandof
   !     J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up*global_aux_var_up%den_kg(iphase)
   !     J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
   !   enddo
      endif
    enddo
  endif
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
                   area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_dn, tor_dn,  dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
        diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
        diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                    weight*rt_parameter%diffusion_coefficient(iphase)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        ! same as dirichlet above
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
          diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
          diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                      weight*rt_parameter%diffusion_coefficient(iphase)
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit
      q = velocity(iphase)
      diffusion = 0D0
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)

      select case(ibndtype)
        case(DIRICHLET_BC)
          if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
            weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
           ! need to account for multiple phases
           ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
            if( iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            if( iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                         weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
           ! same as dirichlet above
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
              diffusion = 0.d0
<<<<<<< local
              if(iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              if(iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                          weight*rt_parameter%diffusion_coefficient(iphase)
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
  endif
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
                             area,rt_parameter,option,velocity,J_dn)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
        diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
        diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                    weight*rt_parameter%diffusion_coefficient(iphase)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
          diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
          diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                      weight*rt_parameter%diffusion_coefficient(iphase)
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit

! super critical CO2 phase
      q = velocity(iphase)
      diffusion = 0D0
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)
  
      select case(ibndtype)
        case(DIRICHLET_BC)
          if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
            weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
           ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
            diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                        weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
            ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
              diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                          weight*rt_parameter%diffusion_coefficient(iphase)
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
        J_dn = J_dn + rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
      else
        J_dn = 0.d0
        do icomp = 1, option%ntrandof
          J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
        enddo
      endif
    enddo
  endif
#endif

end subroutine TBCFluxDerivative
#else
! ************************************************************************** !
!
! TFlux: Computes flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFlux(rt_parameter, &
                 rt_aux_var_up,global_aux_var_up, & 
                 rt_aux_var_dn,global_aux_var_dn, & 
                 coef_up,coef_dn,option,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: coef_up(*), coef_dn(*)
  type(option_type) :: option
  PetscReal :: Res(rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: idof
  PetscInt :: ndof
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: icollcomp
  PetscInt :: icoll
  PetscInt :: iaqcomp

  iphase = 1
  ndof = rt_parameter%naqcomp
  
  ! units = (L water/sec)*(mol/L) = mol/s
  ! total = mol/L water
  Res(1:ndof) = coef_up(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                coef_dn(iphase)*rt_aux_var_dn%total(1:ndof,iphase)

#ifdef REVISED_TRANSPORT
  if (rt_parameter%ncoll > 0) then
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_coll + icoll
      Res(idof) = &
       ! conc_mob = mol/L water
        coef_up(iphase)*rt_aux_var_up%colloid%conc_mob(icoll)+ &
        coef_dn(iphase)*rt_aux_var_dn%colloid%conc_mob(icoll)
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    do icollcomp = 1, rt_parameter%ncollcomp
      iaqcomp = rt_parameter%coll_spec_to_pri_spec(icollcomp)
      ! total_eq_mob = mol/L water
      Res(iaqcomp) = Res(iaqcomp) + &
        coef_up(iphase)*rt_aux_var_up%colloid%total_eq_mob(icollcomp) + &
        coef_dn(iphase)*rt_aux_var_dn%colloid%total_eq_mob(icollcomp)
    enddo
  endif
#endif
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do
     iphase = iphase +1 
     if (iphase > option%nphase) exit
!    super critical CO2 phase have the index 2: need implementation
  
!    units = (L water/sec)*(mol/L) = mol/s
     Res(1:ndof) = Res (1:ndof) + & 
                coef_up(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                coef_dn(iphase)*rt_aux_var_dn%total(1:ndof,iphase)
    enddo
  endif
#endif

end subroutine TFlux

! ************************************************************************** !
!
! TFluxDerivative: Computes derivatives of flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFluxDerivative(rt_parameter, &
                           rt_aux_var_up,global_aux_var_up, & 
                           rt_aux_var_dn,global_aux_var_dn, & 
                           coef_up,coef_dn,option,J_up,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: coef_up(*), coef_dn(*)
  type(option_type) :: option
  PetscReal :: J_up(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_dn(rt_parameter%ncomp,rt_parameter%ncomp)
  
  PetscInt :: iphase
  PetscInt :: icomp
  PetscInt :: icoll
  PetscInt :: idof
  PetscInt :: istart
  PetscInt :: iendaq
 
  iphase = 1
  
  ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  istart = 1
  iendaq = rt_parameter%naqcomp
  if (associated(rt_aux_var_dn%aqueous%dtotal)) then
    J_up(istart:iendaq,istart:iendaq) = rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_up(iphase)
    J_dn(istart:iendaq,istart:iendaq) = rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_dn(iphase)
  else  
    J_up = 0.d0
    J_dn = 0.d0
    do icomp = istart, iendaq
      J_up(icomp,icomp) = coef_up(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_dn(icomp,icomp) = coef_dn(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
    enddo
  endif

#ifdef REVISED_TRANSPORT
  if (rt_parameter%ncoll > 0) then
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_coll + icoll
      J_up(idof,idof) = coef_up(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_dn(idof,idof) = coef_dn(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    ! dRj_dCj - mobile
    ! istart & iend same as above
    J_up(istart:iendaq,istart:iendaq) = J_up(istart:iendaq,istart:iendaq) + &
      rt_aux_var_up%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_up(iphase)
    J_dn(istart:iendaq,istart:iendaq) = J_dn(istart:iendaq,istart:iendaq) + &
      rt_aux_var_dn%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_dn(iphase)
    ! need the below
    ! dRj_dSic
    ! dRic_dSic
    ! dRic_dCj
  endif
#endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit
! super critical CO2 phase

    ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
      if (associated(rt_aux_var_dn%aqueous%dtotal)) then
        J_up(istart:iendaq,istart:iendaq) = J_up(istart:iendaq,istart:iendaq) + &
          rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_up(iphase)
        J_dn(istart:iendaq,istart:iendaq) = J_dn(istart:iendaq,istart:iendaq) + &
          rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_dn(iphase)
      else  
        print *,'Dtotal needed for SC problem. STOP'
        stop 
   !   J_up = 0.d0
   !   J_dn = 0.d0
   !   do icomp = 1, ndof
   !     J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up*global_aux_var_up%den_kg(iphase)
   !     J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
   !   enddo
      endif
    enddo
  endif
#endif

end subroutine TFluxDerivative

#if 0
! ************************************************************************** !
!
! TBCFluxDerivative: Computes derivative of boundary flux term in residual 
!                    function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TBCFluxDerivative(rt_aux_var_up,rt_aux_var_dn, &
                             coef_dn,option,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  PetscReal :: coef_dn(*)
  type(option_type) :: option
  PetscReal :: J_dn(option%ntrandof,option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase

  iphase = 1

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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit

! super critical CO2 phase
      ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
      if (associated(rt_aux_var_dn%dtotal)) then
        J_dn = J_dn + rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
      else
        J_dn = 0.d0
        do icomp = 1, option%ntrandof
          J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
        enddo
      endif
    enddo
  endif
#endif

end subroutine TBCFluxDerivative
#endif
#endif
! ************************************************************************** !
!
! TFluxCoef: Computes blux coefficients for transport matrix
! author: Glenn Hammond
! date: 02/22/10
!
! ************************************************************************** !
subroutine TFluxCoef(option,area,velocity,diffusion,T_up,T_dn)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: velocity(*)
  PetscReal :: diffusion(*)
  PetscReal :: T_up(*), T_dn(*)

  PetscInt :: iphase
  PetscReal :: coef_up, coef_dn
  PetscReal :: q
  
  iphase = 1
  
  q = velocity(iphase)

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion(iphase)+q
    coef_dn = -diffusion(iphase)
  else
    coef_up =  diffusion(iphase)
    coef_dn = -diffusion(iphase)+q
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  T_up(iphase) = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  T_dn(iphase) = coef_dn*area*1000.d0

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do
      iphase = iphase +1 
      if (iphase > option%nphase) exit
    ! super critical CO2 phase have the index 2: need implementation
      q = velocity(iphase)
  
    !upstream weighting
    ! units = (m^3 water/m^2 bulk/sec)
      if (q > 0.d0) then
        coef_up =  diffusion(iphase)+q
        coef_dn = -diffusion(iphase)
      else
        coef_up =  diffusion(iphase)
        coef_dn = -diffusion(iphase)+q
      endif
  
    ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
    !       = L water/sec
      T_up(iphase) = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
      T_dn(iphase) = coef_dn*area*1000.d0
  
    enddo
  endif
#endif

end subroutine TFluxCoef

#ifndef REVISED_TRANSPORT
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
!
! TFluxAdv: Computes advective flux term in residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TFluxAdv(rt_aux_var_up,global_aux_var_up, &
                    rt_aux_var_dn,global_aux_var_dn, &
                    area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  type(reactive_transport_param_type) :: rt_parameter
  PetscReal :: Res(option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: coef_up, coef_dn, sat_up, sat_dn, stp_up, stp_dn, &
               tor_up, tor_dn, por_up, por_dn, dist_up, dist_dn, &
               diffusion, weight
  PetscReal :: q
  
  iphase = 1
  q = velocity(iphase)
 
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up = q
    coef_dn = 0.d0
  else
    coef_up = 0.d0
    coef_dn = q
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do
      iphase = iphase + 1 
      if (iphase > option%nphase) exit
!   super critical CO2 phase have the index 2: need implementation
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

       diffusion = 0.d0 
<<<<<<< local
       if(iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn)+ &
=======
       if(iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn)+ &
>>>>>>> other
                              weight*rt_parameter%diffusion_coefficient(iphase)

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
  endif
#endif

end subroutine TFluxAdv

! ************************************************************************** !
!
! TFluxDerivativeAdv: Computes derivatives of advective flux term in residual 
!                     function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TFluxDerivativeAdv(rt_aux_var_up,global_aux_var_up, &
                              rt_aux_var_dn,global_aux_var_dn, &
                              area,rt_parameter,option,velocity,J_up,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  PetscReal :: J_up(option%ntrandof,option%ntrandof), J_dn(option%ntrandof,option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: q
  PetscReal :: coef_up, coef_dn, sat_up, sat_dn, stp_up, stp_dn, &
               tor_up, tor_dn, por_up, por_dn, dist_up, dist_dn, &
               diffusion, weight
  
  iphase = 1
  q = velocity(iphase)

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up = q
    coef_dn = 0.d0
  else
    coef_up = 0.d0
    coef_dn = q
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
        if(iphase==2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase==2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                                weight*rt_parameter%diffusion_coefficient(iphase)
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
        print *,'Dtotal needed for SC problem. STOP'
        stop 
   !   J_up = 0.d0
   !   J_dn = 0.d0
   !   do icomp = 1, option%ntrandof
   !     J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up*global_aux_var_up%den_kg(iphase)
   !     J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
   !   enddo
      endif
    enddo
  endif
#endif

end subroutine TFluxDerivativeAdv

! ************************************************************************** !
!
! TBCFluxAdv: Computes advective boundary flux term in residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TBCFluxAdv(ibndtype, &
                      rt_aux_var_up,global_aux_var_up, &
                      rt_aux_var_dn,global_aux_var_dn, &
                      area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  PetscReal :: Res(option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: coef_up, coef_dn, sat_up, sat_dn, stp_up, stp_dn, &
               tor_up, tor_dn, por_up, por_dn, dist_up, dist_dn, &
               diffusion, weight
  PetscReal :: q

  iphase = 1
  q = velocity(iphase)
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up = q
    coef_dn = 0.d0
  else
    coef_up = 0.d0
    coef_dn = q
  endif

  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0

  ! units = (L water/sec)*(mol/L) = mol/s  
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)  

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
            if( iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            if( iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                       weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
           ! same as dirichlet above
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
              diffusion = 0.d0
<<<<<<< local
              if(iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn+ &
=======
              if(iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn+ &
>>>>>>> other
                                        weight*rt_parameter%diffusion_coefficient(iphase)
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
  endif
#endif

end subroutine TBCFluxAdv

! ************************************************************************** !
!
! TBCFluxDerivativeAdv: Computes derivative of advective boundary flux term in 
!                       residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TBCFluxDerivativeAdv(ibndtype, &
                                rt_aux_var_up,global_aux_var_up, &
                                rt_aux_var_dn,global_aux_var_dn, &
                                area,rt_parameter,option,velocity,J_dn)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
  type(option_type) :: option
  PetscReal :: J_dn(option%ntrandof,option%ntrandof)
  
  PetscInt :: icomp
  PetscInt :: iphase
  PetscReal :: coef_up, coef_dn, sat_up, sat_dn, stp_up, stp_dn, &
               tor_up, tor_dn, por_up, por_dn, dist_up, dist_dn, &
               diffusion, weight
  PetscReal :: q
  
  iphase = 1
  q = velocity(iphase)
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)  
  if (q > 0.d0) then
    coef_dn = 0.d0
  else
    coef_dn = q
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
            diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                        weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
            ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
              diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                          weight*rt_parameter%diffusion_coefficient(iphase)
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
        J_dn = J_dn + rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
      else
        J_dn = 0.d0
        do icomp = 1, option%ntrandof
          J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
        enddo
      endif
    enddo
  endif
#endif

end subroutine TBCFluxDerivativeAdv

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !

! ************************************************************************** !
!
! TFluxDiff: Computes diffusive flux term in residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TFluxDiff(rt_aux_var_up,global_aux_var_up,por_up,tor_up,dist_up, &
                     rt_aux_var_dn,global_aux_var_dn,por_dn,tor_dn,dist_dn, &
                     area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(option_type) :: option
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
    diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
    diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                weight*rt_parameter%diffusion_coefficient(iphase)
  endif
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion
  endif
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0
  
  ! units = (L water/sec)*(mol/L) = mol/s
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
    do
      iphase = iphase +1 
      if (iphase > option%nphase) exit
! super critical CO2 phase have the index 2: need implementation
      q = velocity(iphase)
      diffusion = 0.d0
      sat_up = global_aux_var_up%sat(iphase)
      sat_dn = global_aux_var_dn%sat(iphase)
  
      if (sat_up > eps .and. sat_dn > eps) then
        stp_up = sat_up*tor_up*por_up 
        stp_dn = sat_dn*tor_dn*por_dn
    ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
        weight = (stp_up*stp_dn)/(stp_up*dist_dn+stp_dn*dist_up)
    ! need to account for multiple phases
    ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec

<<<<<<< local
        if(iphase ==2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase ==2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                              weight*rt_parameter%diffusion_coefficient(iphase)
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
  endif
#endif

end subroutine TFluxDiff

! ************************************************************************** !
!
! TFluxDerivativeDiff: Computes derivatives of diffusive flux term in residual
!                      function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TFluxDerivativeDiff(rt_aux_var_up,global_aux_var_up, &
                               por_up,tor_up,dist_up, &
                               rt_aux_var_dn,global_aux_var_dn, &
                               por_dn,tor_dn,dist_dn, &
                               area,rt_parameter,option,velocity,J_up,J_dn)

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_up, tor_up, dist_up
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
    diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
    diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                weight*rt_parameter%diffusion_coefficient(iphase)
  endif
  
  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
        if(iphase==2) diffusion = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
=======
        if(iphase==2) diffusion = rt_parameter%dispersivity*abs(q)/(dist_up+dist_dn) + &
>>>>>>> other
                                weight*rt_parameter%diffusion_coefficient(iphase)
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
        print *,'Dtotal needed for SC problem. STOP'
        stop 
   !   J_up = 0.d0
   !   J_dn = 0.d0
   !   do icomp = 1, option%ntrandof
   !     J_up(icomp,icomp) = J_up(icomp,icomp) + coef_up*global_aux_var_up%den_kg(iphase)
   !     J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
   !   enddo
      endif
    enddo
  endif
#endif

end subroutine TFluxDerivativeDiff

! ************************************************************************** !
!
! TBCFluxDiff: Computes diffusive boundary flux term in residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TBCFluxDiff(ibndtype, &
                      rt_aux_var_up,global_aux_var_up, &
                      rt_aux_var_dn,global_aux_var_dn, &
                      por_dn,tor_dn,dist_dn, &
                      area,rt_parameter,option,velocity,Res)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  PetscReal :: por_dn, tor_dn,  dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
        diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
        diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                    weight*rt_parameter%diffusion_coefficient(iphase)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        ! same as dirichlet above
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
          diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
          diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                      weight*rt_parameter%diffusion_coefficient(iphase)
        endif    
      endif
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)
  if (q > 0.d0) then
    coef_up =  diffusion
    coef_dn = -diffusion
  else
    coef_up =  diffusion
    coef_dn = -diffusion
  endif

  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  coef_up = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  coef_dn = coef_dn*area*1000.d0

  ! units = (L water/sec)*(mol/L) = mol/s  
  Res(1:option%ntrandof) = coef_up*rt_aux_var_up%total(1:option%ntrandof,iphase) + &
                           coef_dn*rt_aux_var_dn%total(1:option%ntrandof,iphase)  

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
            if( iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            if( iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                         weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
           ! same as dirichlet above
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
              diffusion = 0.d0
<<<<<<< local
              if(iphase == 2) diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              if(iphase == 2) diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                                          weight*rt_parameter%diffusion_coefficient(iphase)
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
  endif
#endif

end subroutine TBCFluxDiff

! ************************************************************************** !
!
! TBCFluxDerivativeDiff: Computes derivative of diffusive boundary flux term 
!                        in residual function
! author: Glenn Hammond
! date: 03/26/09
!
! ************************************************************************** !
subroutine TBCFluxDerivativeDiff(ibndtype, &
                                 rt_aux_var_up,global_aux_var_up, &
                                 rt_aux_var_dn,global_aux_var_dn, &
                                 por_dn,tor_dn,dist_dn, &
                                 area,rt_parameter,option,velocity,J_dn)

  use Option_module

  implicit none
  
  PetscInt :: ibndtype
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: por_dn, tor_dn, dist_dn
  PetscReal :: area
  PetscReal :: velocity(1)
  type(reactive_transport_param_type) :: rt_parameter
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
<<<<<<< local
        diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
        diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                    weight*rt_parameter%diffusion_coefficient(iphase)
      endif    
    case(DIRICHLET_ZERO_GRADIENT_BC)
      if (q >= 0.d0) then
        if (sat_up > eps .and. sat_dn > eps) then
          ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
          weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
          ! need to account for multiple phases
          ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
          diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
          diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                      weight*rt_parameter%diffusion_coefficient(iphase)
        endif  
      endif  
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select

  !upstream weighting
  ! units = (m^3 water/m^2 bulk/sec)  
  if (q > 0.d0) then
    coef_dn = -diffusion
  else
    coef_dn = -diffusion
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
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE) then
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
<<<<<<< local
            diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
            diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                        weight*rt_parameter%diffusion_coefficient(iphase)
          endif    
        case(DIRICHLET_ZERO_GRADIENT_BC)
          if (q >= 0.d0) then
            if (sat_up > eps .and. sat_dn > eps) then
            ! units = (m^3 water/m^3 por)*(m^3 por/m^3 bulk)/(m bulk) = m^3 water/m^4 bulk 
              weight = tor_dn*por_dn*(sat_up*sat_dn)/((sat_up+sat_dn)*dist_dn)
            ! need to account for multiple phases
            ! units = (m^3 water/m^4 bulk)*(m^2 bulk/sec) = m^3 water/m^2 bulk/sec
<<<<<<< local
              diffusion = rt_parameter%dispersivity*dabs(q)/dist_dn + &
=======
              diffusion = rt_parameter%dispersivity*abs(q)/dist_dn + &
>>>>>>> other
                          weight*rt_parameter%diffusion_coefficient(iphase)
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
        J_dn = J_dn + rt_aux_var_dn%dtotal(:,:,iphase)*coef_dn*1000.d0
      else
        J_dn = 0.d0
        do icomp = 1, option%ntrandof
          J_dn(icomp,icomp) = J_dn(icomp,icomp) + coef_dn*global_aux_var_dn%den_kg(iphase)
        enddo
      endif
    enddo
  endif
#endif

end subroutine TBCFluxDerivativeDiff
#endif

end module Transport_module
