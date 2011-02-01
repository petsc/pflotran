module Transport_module

  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Matrix_Block_Aux_module  

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
            TFluxCoef, &
            TSrcSinkCoef, &
            TFlux_CD, &
            TFluxDerivative_CD, &
            TFluxCoef_CD
              
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
    diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
                        weight*rt_parameter%diffusion_coefficient(iphase)
  endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
       .or. option%iflowmode == FLASH2_MODE) then
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
        if(iphase ==2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/(dist_up+dist_dn) + &
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
        diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
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
          diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
                              weight*rt_parameter%diffusion_coefficient(iphase)
        endif    
      endif
    case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
  end select



! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
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
            if( iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
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
              if(iphase == 2) diffusion(iphase) = rt_parameter%dispersivity*dabs(q)/dist_dn + &
                                        weight*rt_parameter%diffusion_coefficient(iphase)
            endif    
          endif
        case(CONCENTRATION_SS,NEUMANN_BC,ZERO_GRADIENT_BC)
      end select
    enddo
  endif
#endif

end subroutine TDiffusionBC


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
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
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
! TFlux: Computes flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFlux_CD(rt_parameter, &
                 rt_aux_var_up,global_aux_var_up, & 
                 rt_aux_var_dn,global_aux_var_dn, & 
                 coef_11,coef_12,coef_21,coef_22,option,Res_1,Res_2)

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: coef_11(*), coef_12(*), coef_21(*), coef_22(*)
  type(option_type) :: option
  PetscReal :: Res_1(rt_parameter%ncomp)
  PetscReal :: Res_2(rt_parameter%ncomp)
  
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
  Res_1(1:ndof) = coef_11(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                  coef_12(iphase)*rt_aux_var_dn%total(1:ndof,iphase)
  Res_2(1:ndof) = coef_21(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                  coef_22(iphase)*rt_aux_var_dn%total(1:ndof,iphase)

  if (rt_parameter%ncoll > 0) then
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_coll + icoll
       ! conc_mob = mol/L water
      Res_1(idof) = coef_11(iphase)*rt_aux_var_up%colloid%conc_mob(icoll)+ &
                    coef_12(iphase)*rt_aux_var_dn%colloid%conc_mob(icoll)
      Res_2(idof) = coef_21(iphase)*rt_aux_var_up%colloid%conc_mob(icoll)+ &
                    coef_22(iphase)*rt_aux_var_dn%colloid%conc_mob(icoll)
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    do icollcomp = 1, rt_parameter%ncollcomp
      iaqcomp = rt_parameter%coll_spec_to_pri_spec(icollcomp)
      ! total_eq_mob = mol/L water
      Res_1(iaqcomp) = Res_1(iaqcomp) + &
        coef_11(iphase)*rt_aux_var_up%colloid%total_eq_mob(icollcomp) + &
        coef_12(iphase)*rt_aux_var_dn%colloid%total_eq_mob(icollcomp)
      Res_2(iaqcomp) = Res_2(iaqcomp) + &
        coef_21(iphase)*rt_aux_var_up%colloid%total_eq_mob(icollcomp) + &
        coef_22(iphase)*rt_aux_var_dn%colloid%total_eq_mob(icollcomp)
    enddo
  endif
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
    do
     iphase = iphase +1 
     if (iphase > option%nphase) exit
!    super critical CO2 phase have the index 2: need implementation
  
!    units = (L water/sec)*(mol/L) = mol/s
     Res_1(1:ndof) = Res_1(1:ndof) + &
                       coef_11(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                       coef_12(iphase)*rt_aux_var_dn%total(1:ndof,iphase)
     Res_2(1:ndof) = Res_2(1:ndof) + &
                       coef_21(iphase)*rt_aux_var_up%total(1:ndof,iphase) + &
                       coef_22(iphase)*rt_aux_var_dn%total(1:ndof,iphase)
    enddo
  endif
#endif

end subroutine TFlux_CD

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

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
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

! ************************************************************************** !
!
! TFluxDerivative: Computes derivatives of flux term in residual function
! author: Glenn Hammond
! date: 02/15/08
!
! ************************************************************************** !
subroutine TFluxDerivative_CD(rt_parameter, &
                           rt_aux_var_up,global_aux_var_up, & 
                           rt_aux_var_dn,global_aux_var_dn, & 
                           coef_11,coef_12,coef_21,coef_22,option, &
                           J_11,J_12,J_21,J_22)

  use Option_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_aux_var_up, rt_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn 
  PetscReal :: coef_11(*), coef_12(*), coef_21(*), coef_22(*)
  type(option_type) :: option
  PetscReal :: J_11(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_12(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_21(rt_parameter%ncomp,rt_parameter%ncomp), &
               J_22(rt_parameter%ncomp,rt_parameter%ncomp)
  
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
    J_11(istart:iendaq,istart:iendaq) = rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_11(iphase)
    J_12(istart:iendaq,istart:iendaq) = rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_12(iphase)
    J_21(istart:iendaq,istart:iendaq) = rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_21(iphase)
    J_22(istart:iendaq,istart:iendaq) = rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_22(iphase)
  else  
    J_11 = 0.d0
    J_12 = 0.d0
    J_21 = 0.d0
    J_22 = 0.d0
    do icomp = istart, iendaq
      J_11(icomp,icomp) = coef_11(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_12(icomp,icomp) = coef_12(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
      J_21(icomp,icomp) = coef_21(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_22(icomp,icomp) = coef_22(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
    enddo
  endif

  if (rt_parameter%ncoll > 0) then
    do icoll = 1, rt_parameter%ncoll
      idof = rt_parameter%offset_coll + icoll
      J_11(idof,idof) = coef_11(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_12(idof,idof) = coef_12(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
      J_21(idof,idof) = coef_21(iphase)*global_aux_var_up%den_kg(iphase)*1.d-3
      J_22(idof,idof) = coef_22(iphase)*global_aux_var_dn%den_kg(iphase)*1.d-3
    enddo
  endif
  if (rt_parameter%ncollcomp > 0) then
    ! dRj_dCj - mobile
    ! istart & iend same as above
    J_11(istart:iendaq,istart:iendaq) = J_11(istart:iendaq,istart:iendaq) + &
      rt_aux_var_up%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_11(iphase)
    J_12(istart:iendaq,istart:iendaq) = J_12(istart:iendaq,istart:iendaq) + &
      rt_aux_var_dn%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_12(iphase)
    J_21(istart:iendaq,istart:iendaq) = J_21(istart:iendaq,istart:iendaq) + &
      rt_aux_var_up%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_21(iphase)
    J_22(istart:iendaq,istart:iendaq) = J_22(istart:iendaq,istart:iendaq) + &
      rt_aux_var_dn%colloid%dRj_dCj%dtotal(:,:,iphase)*coef_22(iphase)
    ! need the below
    ! dRj_dSic
    ! dRic_dSic
    ! dRic_dCj
  endif

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
    do 
      iphase = iphase + 1
      if (iphase > option%nphase) exit
! super critical CO2 phase

    ! units = (m^3 water/sec)*(kg water/L water)*(1000L water/m^3 water) = kg water/sec
      if (associated(rt_aux_var_dn%aqueous%dtotal)) then
        J_11(istart:iendaq,istart:iendaq) = J_11(istart:iendaq,istart:iendaq) + &
          rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_11(iphase)
        J_12(istart:iendaq,istart:iendaq) = J_12(istart:iendaq,istart:iendaq) + &
          rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_12(iphase)
        J_21(istart:iendaq,istart:iendaq) = J_21(istart:iendaq,istart:iendaq) + &
          rt_aux_var_up%aqueous%dtotal(:,:,iphase)*coef_21(iphase)
        J_22(istart:iendaq,istart:iendaq) = J_22(istart:iendaq,istart:iendaq) + &
          rt_aux_var_dn%aqueous%dtotal(:,:,iphase)*coef_22(iphase)
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

end subroutine TFluxDerivative_CD

! ************************************************************************** !
!
! TFluxCoef: Computes flux coefficients for transport matrix
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

  if (option%use_upwinding) then
    ! upstream weighting
    ! units = (m^3 water/m^2 bulk/sec)
    if (q > 0.d0) then
      coef_up =  diffusion(iphase)+q
      coef_dn = -diffusion(iphase)
    else
      coef_up =  diffusion(iphase)
      coef_dn = -diffusion(iphase)+q
    endif
  else
    ! central difference, currently assuming uniform grid spacing
    ! units = (m^3 water/m^2 bulk/sec)
    ! 
    ! coef_up/coef_dn needs to flip sign of dispersion, but preserve sign of
    ! of advection.  Since I simply flip the sign of coef_up/coef_dn outside
    ! in RTResidual and RTJacobian, flipping dispersion without flipping
    ! advection is not possible, and thus central difference will not work
    ! for now.
    if (q > 0.d0) then
      coef_up =  diffusion(iphase)+ 0.5d0*q
      coef_dn = -diffusion(iphase)+ 0.5d0*q
      ! coef_up_dif = diffusion(iphase)
      ! coef_up_adv = -0.5d0*q
      ! coef_dn_dif = -diffusion(iphase)
      ! coef_dn_adv = 0.5d0*q
    else
      coef_up =  diffusion(iphase)+ 0.5d0*q
      coef_dn = -diffusion(iphase)+ 0.5d0*q
    endif
  endif  
  
  ! units = (m^3 water/m^2 bulk/sec)*(m^2 bulk)*(1000 L water/m^3 water)
  !       = L water/sec
  T_up(iphase) = coef_up*area*1000.d0  ! 1000 converts m^3 -> L
  T_dn(iphase) = coef_dn*area*1000.d0

! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
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

! ************************************************************************** !
!
! TFluxCoef_CD: Computes flux coefficients for transport matrix
! author: Glenn Hammond
! date: 02/22/10
!
! ************************************************************************** !
subroutine TFluxCoef_CD(option,area,velocity,diffusion,T_11,T_12,T_21,T_22)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: velocity(*)
  PetscReal :: diffusion(*)
  PetscReal :: T_11(*), T_12(*), T_21(*), T_22(*)

  PetscInt :: iphase
  PetscReal :: coef_up, coef_dn
  PetscReal :: tempreal
  PetscReal :: weight
  PetscReal :: q
  
  iphase = 1
  
  ! T_11 = diagonal term for upwind cell (row)
  ! T_12 = off diagonal term for upwind cell (row)
  ! T_21 = off diagonal term for downwind cell (row)
  ! T_22 = diagonal term for downwind cell (row)
  
  q = velocity(iphase)

  if (option%use_upwinding) then
    ! upstream weighting
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
    tempreal = area*1000.d0
    T_11(iphase) = coef_up*tempreal  ! 1000 converts m^3 -> L
    T_12(iphase) = coef_dn*tempreal
    T_21(iphase) = -T_11(iphase)
    T_22(iphase) = -T_12(iphase)

  else
    ! central difference, currently assuming uniform grid spacing
    ! units = (m^3 water/m^2 bulk/sec)
    ! 
    weight = 0.5d0
    tempreal = area*1000.d0
    T_11(iphase) = (diffusion(iphase) + weight*q)*tempreal
    T_12(iphase) = (-diffusion(iphase) + weight*q)*tempreal
    T_21(iphase) = -(diffusion(iphase) + weight*q)*tempreal
    T_22(iphase) = (diffusion(iphase) - weight*q)*tempreal
  endif  
  
! Add in multiphase, clu 12/29/08
#ifdef CHUAN_CO2  
  if (option%iflowmode == MPH_MODE .or. option%iflowmode == IMS_MODE &
      .or. option%iflowmode == FLASH2_MODE) then
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

end subroutine TFluxCoef_CD

! ************************************************************************** !
!
! TSrcSinkCoef: Computes src/sink coefficients for transport matrix
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine TSrcSinkCoef(option,qsrc,flow_src_sink_type,tran_src_sink_type, &
                        por,sat,vol,den,scale,kg_per_sec,T_in,T_out)

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: qsrc
  PetscInt :: flow_src_sink_type
  PetscInt :: tran_src_sink_type
  PetscReal :: por
  PetscReal :: sat
  PetscReal :: vol
  PetscReal :: den
  PetscReal :: scale
  PetscBool :: kg_per_sec
  PetscReal :: T_in ! coefficient that scales concentration at cell
  PetscReal :: T_out ! concentration that scales external concentration
      
  PetscReal :: rate
  
  T_in = 0.d0 
  T_out = 0.d0
     
  select case(tran_src_sink_type)
    case(EQUILIBRIUM_SS)
      ! units should be mol/sec
      rate = 1.d-6 ! units 1/sec
      T_in = rate*por*sat*vol*1000.d0 ! units L water/sec
      T_out = -1.d0*T_in
    case(MASS_RATE_SS)
      ! in this case, rt_auxvar_bc%total actually holds the mass rate
      T_in = 0.d0
      T_out = -1.d0
    case default
      if (qsrc > 0.d0) then ! injection
        T_in = 0.d0
        select case(flow_src_sink_type)
          case(MASS_RATE_SS)
            T_out = -1.d0*qsrc/den*1000.d0 ! kg water/sec / kg water/m^3 * 1000 L/m^3 -> L/sec
          case(SCALED_MASS_RATE_SS)
            T_out = -1.d0*qsrc/den*1000.d0*scale ! m^3/sec * 1000 m^3/L -> L/s
          case(VOLUMETRIC_RATE_SS)
            T_out = -1.d0*qsrc*1000.d0 ! m^3/sec * 1000 m^3/L -> L/s
          case(SCALED_VOLUMETRIC_RATE_SS)
            T_out = -1.d0*qsrc*1000.d0*scale ! m^3/sec * 1000 m^3/L -> L/s
        end select
      else
        T_out = 0.d0
        select case(flow_src_sink_type)
          case(MASS_RATE_SS)
            T_in = -1.d0*qsrc/den*1000.d0 ! kg water/sec / kg water/m^3 * 1000 L/m^3 -> L/sec
          case(SCALED_MASS_RATE_SS)
            T_in = -1.d0*qsrc/den*1000.d0*scale ! m^3/sec * 1000 m^3/L -> L/s
          case(VOLUMETRIC_RATE_SS)
            T_in = -1.d0*qsrc*1000.d0 ! m^3/sec * 1000 m^3/L -> L/s
          case(SCALED_VOLUMETRIC_RATE_SS)
            T_in = -1.d0*qsrc*1000.d0*scale ! m^3/sec * 1000 m^3/L -> L/s
        end select
      endif
  end select

  ! up to this point, units = L water/sec, which is correct for residual function
  ! for Jacobian, need kg water/sec
  if (kg_per_sec) then
    T_in = T_in * den / 1000.d0
    T_out = T_out * den / 1000.d0
  endif
  
end subroutine TSrcSinkCoef
  
end module Transport_module
