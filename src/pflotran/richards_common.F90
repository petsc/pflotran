module Richards_Common_module

  use Richards_Aux_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6
  
  public RichardsAccumulation, &
         RichardsAccumDerivative, &
         RichardsFlux, &
         RichardsFluxDerivative, &
         RichardsBCFlux, &
         RichardsBCFluxDerivative
  
contains

! ************************************************************************** !
!
! RichardsAccumDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsAccumDerivative(rich_aux_var,global_aux_var,por,vol, &
                                   option,sat_func,J)

  use Option_module
  use Saturation_Function_module
  
  implicit none

  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: vol, por
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: vol_over_dt
  PetscReal :: tempreal 

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert
  type(global_auxvar_type) :: global_aux_var_pert
  PetscReal :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  vol_over_dt = vol/option%flow_dt
      
!#define USE_COMPRESSIBLITY
#ifndef USE_COMPRESSIBLITY  
  ! accumulation term units = dkmol/dp
  J(1,1) = (global_aux_var%sat(1)*rich_aux_var%dden_dp+ &
            rich_aux_var%dsat_dp*global_aux_var%den(1))* &
           por*vol_over_dt

#else
  tempreal = exp(-1.d-7*(global_aux_var%pres(1)-option%reference_pressure))
  J(1,1) = ((global_aux_var%sat(1)*rich_aux_var%dden_dp+ &
             rich_aux_var%dsat_dp*global_aux_var%den(1))* &
            1.d0-(1.d0-por)*tempreal + &
            global_aux_var%sat(1)*global_aux_var%den(1)* &
            (por-1.d0)*-1.d-7*tempreal)* &
           vol_over_dt
#endif  
  
  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert,option)  
    call RichardsAuxVarCopy(rich_aux_var,rich_aux_var_pert,option)
    call GlobalAuxVarCopy(global_aux_var,global_aux_var_pert,option)
    x(1) = global_aux_var%pres(1)
    call RichardsAccumulation(rich_aux_var,global_aux_var,por,vol,option,res)
    ideriv = 1
    pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
    x_pert = x
    if (x_pert(ideriv) < option%reference_pressure) pert = -1.d0*pert
    x_pert(ideriv) = x_pert(ideriv) + pert
    
    call RichardsAuxVarCompute(x_pert(1),rich_aux_var_pert,global_aux_var_pert, &
                               sat_func,0.d0,0.d0,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', aux_var%dvis_dp, (aux_var_pert%vis-aux_var%vis)/pert(ideriv)
!         print *, 'dkr_dp:', aux_var%dkr_dp, (aux_var_pert%kr-aux_var%kr)/pert(ideriv)
          print *, 'dsat_dp:', aux_var%dsat_dp, (global_aux_var_pert%sat-global_aux_var%sat)/pert
          print *, 'dden_dp:', aux_var%dden_dp, (global_aux_var_pert%den-global_aux_var%den)/pert
!          print *, 'dkvr_dp:', aux_var%dkvr_dp, (rich_aux_var_pert%kvr-rich_aux_var%kvr)/pert
          print *, 'dkvr_x_dp:', aux_var%dkvr_x_dp, (rich_aux_var_pert%kvr_x-rich_aux_var%kvr_x)/pert
          print *, 'dkvr_y_dp:', aux_var%dkvr_y_dp, (rich_aux_var_pert%kvr_y-rich_aux_var%kvr_y)/pert
          print *, 'dkvr_z_dp:', aux_var%dkvr_z_dp, (rich_aux_var_pert%kvr_z-rich_aux_var%kvr_z)/pert
      end select     
#endif     
    call RichardsAccumulation(rich_aux_var_pert,global_aux_var_pert,por,vol, &
                              option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
    call GlobalAuxVarStrip(global_aux_var_pert)  
  endif
   
end subroutine RichardsAccumDerivative

! ************************************************************************** !
!
! RichardsAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsAccumulation(rich_aux_var,global_aux_var,por,vol, &
                                option,Res)

  use Option_module
  
  implicit none

  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: vol, por, por1
       
  ! accumulation term units = kmol/s
#ifndef USE_COMPRESSIBILITY
    Res(1) = global_aux_var%sat(1) * global_aux_var%den(1) * por * vol / &
           option%flow_dt
#else
    por1 = 1.d0-(1.d0-por)*exp(-1.d-7*(global_aux_var%pres(1)- &
                                       option%reference_pressure))
    Res(1) = global_aux_var%sat(1) * global_aux_var%den(1) * por1 * vol / &
           option%flow_dt
#endif
    
end subroutine RichardsAccumulation

! ************************************************************************** !
!
! RichardsFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsFluxDerivative(rich_aux_var_up,global_aux_var_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_aux_var_dn,global_aux_var_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  option,sat_func_up,sat_func_dn,Jup,Jdn)
  use Option_module 
  use Saturation_Function_module                        
  
  implicit none
  
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy, area, dist(3)
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
  PetscReal :: q
  PetscReal :: ukvr,Dq
!  PetscReal :: ukvr_x, ukvr_y, ukvr_z, Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
!  PetscReal :: dukvr_x_dp_up, dukvr_x_dp_dn
!  PetscReal :: dukvr_y_dp_up, dukvr_y_dp_dn
!  PetscReal :: dukvr_z_dp_up, dukvr_z_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert_up, rich_aux_var_pert_dn
  type(global_auxvar_type) :: global_aux_var_pert_up, global_aux_var_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
  v_darcy = 0.D0 
  ukvr = 0.d0
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)
    dden_ave_dp_up = upweight*rich_aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_aux_var_dn%dden_dp

    gravity = (upweight*global_aux_var_up%den(1) + &
              (1.D0-upweight)*global_aux_var_dn%den(1)) &
              * FMWH2O * dist_gravity
    dgravity_dden_up = upweight*FMWH2O*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*rich_aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_aux_var_dn%dden_dp

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dist(1))==1) then
        ukvr = rich_aux_var_up%kvr_x
        dukvr_dp_up = rich_aux_var_up%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_aux_var_up%kvr_y
        dukvr_dp_up = rich_aux_var_up%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_aux_var_up%kvr_z
        dukvr_dp_up = rich_aux_var_up%dkvr_z_dp
      end if
#else
      ukvr = rich_aux_var_up%kvr
      dukvr_dp_up = rich_aux_var_up%dkvr_dp
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY    
      if (dabs(dist(1))==1) then
        ukvr = rich_aux_var_dn%kvr_x
        dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
      else if (dabs(dist(2))==1) then
        ukvr = rich_aux_var_dn%kvr_y
        dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
      else if (dabs(dist(3))==1) then
        ukvr = rich_aux_var_dn%kvr_z
        dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
      end if
#else
      ukvr = rich_aux_var_dn%kvr
      dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
    endif
  endif

 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call RichardsAuxVarCopy(rich_aux_var_up,rich_aux_var_pert_up,option)
    call RichardsAuxVarCopy(rich_aux_var_dn,rich_aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    x_up(1) = global_aux_var_up%pres(1)
    x_dn(1) = global_aux_var_dn%pres(1)
    call RichardsFlux(rich_aux_var_up,global_aux_var_up,por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_dn,global_aux_var_dn,por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res)
    ideriv = 1
!    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_up = max(dabs(x_up(ideriv)*perturbation_tolerance),0.1d0)
    if (x_up(ideriv) < option%reference_pressure) pert_up = -1.d0*pert_up
!    pert_dn = x_dn(ideriv)*perturbation_tolerance
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_up = x_up
    x_pert_dn = x_dn
    x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    call RichardsAuxVarCompute(x_pert_up(1),rich_aux_var_pert_up, &
                               global_aux_var_pert_up,sat_func_up, &
                               0.d0,0.d0,option)
    call RichardsAuxVarCompute(x_pert_dn(1),rich_aux_var_pert_dn, &
                               global_aux_var_pert_dn,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsFlux(rich_aux_var_pert_up,global_aux_var_pert_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_dn,global_aux_var_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res_pert_up)
    call RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_aux_var_pert_dn,global_aux_var_pert_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)    
  endif

end subroutine RichardsFluxDerivative

! ************************************************************************** !
!
! RichardsFlux: Computes the internal flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsFlux(rich_aux_var_up,global_aux_var_up, &
                        por_up,sir_up,dd_up,perm_up, &
                        rich_aux_var_dn,global_aux_var_dn, &
                        por_dn,sir_dn,dd_dn,perm_dn, &
                        area, dist, dist_gravity,upweight, &
                        option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy,area, dist(3)
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
  fluxm = 0.d0
  v_darcy = 0.D0  
  ukvr = 0.d0
  
! Flow term
  if (global_aux_var_up%sat(1) > sir_up .or. global_aux_var_dn%sat(1) > sir_dn) then
    if (global_aux_var_up%sat(1) <eps) then 
      upweight=0.d0
    else if (global_aux_var_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    

    density_ave = upweight*global_aux_var_up%den(1)+ &
                  (1.D0-upweight)*global_aux_var_dn%den(1)

    gravity = (upweight*global_aux_var_up%den(1) + &
              (1.D0-upweight)*global_aux_var_dn%den(1)) &
              * FMWH2O * dist_gravity

    dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1)  + gravity

    if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_up%kvr_z
      end if
#else
      ukvr = rich_aux_var_up%kvr
#endif
    else
#ifdef USE_ANISOTROPIC_MOBILITY       
      if (dabs(dabs(norm(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(norm(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(norm(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
#else
      ukvr = rich_aux_var_dn%kvr
#endif
    endif      

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area

      fluxm = q*density_ave       
    endif
  endif 

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsFlux

! ************************************************************************** !
!
! RichardsBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsBCFluxDerivative(ibndtype,aux_vars, &
                                    rich_aux_var_up,global_aux_var_up, &
                                    rich_aux_var_dn,global_aux_var_dn, &
                                    por_dn,sir_dn,perm_dn, &
                                    area, dist,option, &
                                    sat_func_dn,Jdn)
  use Option_module
  use Saturation_Function_module
#ifdef SURFACE_FLOW
  use Water_EOS_module
#endif
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn
  PetscReal :: area
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(0:3)
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscReal :: v_darcy
  PetscReal :: q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: dden_ave_dp_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn
  PetscReal :: dukvr_dp_dn
  PetscReal :: dq_dp_dn
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  type(richards_auxvar_type) :: rich_aux_var_pert_dn, rich_aux_var_pert_up
  type(global_auxvar_type) :: global_aux_var_pert_dn, global_aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
#ifdef SURFACE_FLOW
  PetscReal :: rho, v_darcy_allowable
#endif
  
  v_darcy = 0.d0
  ukvr = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_dn = 0.d0
        
  ! Flow
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC,HET_SURF_SEEPAGE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = aux_vars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_aux_var_up%den(1) + (1.D0-upweight)*global_aux_var_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*rich_aux_var_dn%dden_dp

        gravity = (upweight*global_aux_var_up%den(1) + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)) &
                  * FMWH2O * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*rich_aux_var_dn%dden_dp

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC .or. &
            pressure_bc_type == HET_SURF_SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY  
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_x
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_y
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_aux_var_up%kvr_z
          end if
#else
          ukvr = rich_aux_var_up%kvr
#endif
        else
#ifdef USE_ANISOTROPIC_MOBILITY
          if (dabs(dabs(dist(1))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_x
            dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
          else if (dabs(dabs(dist(2))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_y
            dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
          else if (dabs(dabs(dist(3))-1) < 1e-6) then
            ukvr = rich_aux_var_dn%kvr_z
            dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
          end if
#else
          ukvr = rich_aux_var_dn%kvr
          dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
#ifdef SURFACE_FLOW
          ! If running with surface-flow model, ensure (darcy_velocity*dt) does
          ! not exceed depth of standing water.
          if(pressure_bc_type == HET_SURF_SEEPAGE_BC .and. option%nsurfflowdof>0) then
            call density(option%reference_temperature,option%reference_pressure,rho)
            v_darcy_allowable = (global_aux_var_up%pres(1)-option%reference_pressure) &
                                /option%flow_dt/(-option%gravity(3))/rho
            if(v_darcy > v_darcy_allowable) then
              dphi_dp_dn = 0.d0
              v_darcy = v_darcy_allowable
            endif
          endif
#endif
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi + ukvr*dphi_dp_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
        else 
          density_ave = global_aux_var_dn%den(1)
          dden_ave_dp_dn = rich_aux_var_dn%dden_dp
        endif 
        q = v_darcy * area
      endif

    case(UNIT_GRADIENT_BC)

      if (global_aux_var_dn%sat(1) > sir_dn) then

        dphi = dot_product(option%gravity,dist(1:3))*global_aux_var_dn%den(1)*FMWH2O      
        density_ave = global_aux_var_dn%den(1)

        dphi_dp_dn = dot_product(option%gravity,dist(1:3))*rich_aux_var_dn%dden_dp*FMWH2O
        dden_ave_dp_dn = rich_aux_var_dn%dden_dp

        ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
        if (dabs(dabs(dist(1))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_x
          dukvr_dp_dn = rich_aux_var_dn%dkvr_x_dp
        else if (dabs(dabs(dist(2))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_y
          dukvr_dp_dn = rich_aux_var_dn%dkvr_y_dp
        else if (dabs(dabs(dist(3))-1) < 1e-6) then
          ukvr = rich_aux_var_dn%kvr_z
          dukvr_dp_dn = rich_aux_var_dn%dkvr_z_dp
        end if
#else
        ukvr = rich_aux_var_dn%kvr
        dukvr_dp_dn = rich_aux_var_dn%dkvr_dp
#endif
     
        if (ukvr*perm_dn>floweps) then
          v_darcy = perm_dn * ukvr * dphi
          q = v_darcy*area
          dq_dp_dn = perm_dn*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif 
     
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  if (option%numerical_derivatives_flow) then
    call GlobalAuxVarInit(global_aux_var_pert_up,option)
    call GlobalAuxVarInit(global_aux_var_pert_dn,option)  
    call RichardsAuxVarCopy(rich_aux_var_up,rich_aux_var_pert_up,option)
    call RichardsAuxVarCopy(rich_aux_var_dn,rich_aux_var_pert_dn,option)
    call GlobalAuxVarCopy(global_aux_var_up,global_aux_var_pert_up,option)
    call GlobalAuxVarCopy(global_aux_var_dn,global_aux_var_pert_dn,option)
    x_up(1) = global_aux_var_up%pres(1)
    x_dn(1) = global_aux_var_dn%pres(1)
    ideriv = 1
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_up(ideriv) = x_dn(ideriv)
    endif
    call RichardsBCFlux(ibndtype,aux_vars, &
                        rich_aux_var_up,global_aux_var_up, &
                        rich_aux_var_dn,global_aux_var_dn, &
                        por_dn,sir_dn,perm_dn, &
                        area,dist,option,v_darcy,res)
    if (pressure_bc_type == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    ideriv = 1
!    pert_dn = x_dn(ideriv)*perturbation_tolerance    
    pert_dn = max(dabs(x_dn(ideriv)*perturbation_tolerance),0.1d0)
    if (x_dn(ideriv) < option%reference_pressure) pert_dn = -1.d0*pert_dn
    x_pert_dn = x_dn
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    x_pert_up = x_up
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_pert_up(ideriv) = x_pert_dn(ideriv)
    endif   
    call RichardsAuxVarCompute(x_pert_dn(1),rich_aux_var_pert_dn, &
                               global_aux_var_pert_dn,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsAuxVarCompute(x_pert_up(1),rich_aux_var_pert_up, &
                               global_aux_var_pert_up,sat_func_dn, &
                               0.d0,0.d0,option)
    call RichardsBCFlux(ibndtype,aux_vars, &
                        rich_aux_var_pert_up,global_aux_var_pert_up, &
                        rich_aux_var_pert_dn,global_aux_var_pert_dn, &
                        por_dn,sir_dn,perm_dn, &
                        area,dist,option,v_darcy,res_pert_dn)
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_aux_var_pert_up)
    call GlobalAuxVarStrip(global_aux_var_pert_dn)      
  endif

end subroutine RichardsBCFluxDerivative

! ************************************************************************** !
!
! RichardsBCFlux: Computes the  boundary flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsBCFlux(ibndtype,aux_vars, &
                          rich_aux_var_up, global_aux_var_up, &
                          rich_aux_var_dn, global_aux_var_dn, &
                          por_dn, sir_dn, perm_dn, &
                          area, dist, option,v_darcy,Res)
  use Option_module
#ifdef SURFACE_FLOW
  use Water_EOS_module
#endif
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(richards_auxvar_type) :: rich_aux_var_up, rich_aux_var_dn
  type(global_auxvar_type) :: global_aux_var_up, global_aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn
  PetscReal :: v_darcy, area
  ! dist(0) = magnitude
  ! dist(1:3) = unit vector
  ! dist(0)*dist(1:3) = vector
  PetscReal :: dist(0:3)
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm,q,density_ave
  PetscReal :: ukvr,diffdp,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscInt :: pressure_bc_type
  PetscReal :: dphi_x,dphi_y,dphi_z
#ifdef SURFACE_FLOW
  PetscReal :: rho, v_darcy_allowable
#endif
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  ukvr = 0.d0

  ! Flow  
  pressure_bc_type = ibndtype(RICHARDS_PRESSURE_DOF)
  select case(pressure_bc_type)
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC,HET_SURF_SEEPAGE_BC)

      ! dist(0) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3) = vector(3) - unit vector
      dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

      if (pressure_bc_type == CONDUCTANCE_BC) then
        Dq = aux_vars(RICHARDS_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dist(0)
      endif
      ! Flow term
      if (global_aux_var_up%sat(1) > sir_dn .or. global_aux_var_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_aux_var_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_aux_var_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_aux_var_up%den(1)+(1.D0-upweight)*global_aux_var_dn%den(1)
   
        gravity = (upweight*global_aux_var_up%den(1) + &
                  (1.D0-upweight)*global_aux_var_dn%den(1)) &
                  * FMWH2O * dist_gravity
       
        dphi = global_aux_var_up%pres(1) - global_aux_var_dn%pres(1) + gravity

        if (pressure_bc_type == SEEPAGE_BC .or. &
            pressure_bc_type == CONDUCTANCE_BC .or. &
            pressure_bc_type == HET_SURF_SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. global_aux_var_up%pres(1)-option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
   
       if (dphi>=0.D0) then
#ifdef USE_ANISOTROPIC_MOBILITY       
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_up%kvr_z
         end if
#else
         ukvr = rich_aux_var_up%kvr
#endif
       else
#ifdef USE_ANISOTROPIC_MOBILITY
         if (dabs(dabs(dist(1))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_x
         else if (dabs(dabs(dist(2))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_y
         else if (dabs(dabs(dist(3))-1) < 1e-6) then
           ukvr = rich_aux_var_dn%kvr_z
         end if
#else
         ukvr = rich_aux_var_dn%kvr
#endif
       endif
        
       if (ukvr*Dq>floweps) then
        v_darcy = Dq * ukvr * dphi
#ifdef SURFACE_FLOW
          ! If running with surface-flow model, ensure (darcy_velocity*dt) does
          ! not exceed depth of standing water.
          if(pressure_bc_type == HET_SURF_SEEPAGE_BC .and. option%nsurfflowdof>0) then
            call density(option%reference_temperature,option%reference_pressure,rho)
            v_darcy_allowable = (global_aux_var_up%pres(1)-option%reference_pressure) &
                                /option%flow_dt/(-option%gravity(3))/rho
            if (v_darcy > v_darcy_allowable) v_darcy = v_darcy_allowable
          endif
#endif
       endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)

        if (v_darcy > 0.d0) then 
          density_ave = global_aux_var_up%den(1)
        else 
          density_ave = global_aux_var_dn%den(1)
        endif 
      endif

    case(UNIT_GRADIENT_BC)

      dphi = dot_product(option%gravity,dist(1:3))*global_aux_var_dn%den(1)*FMWH2O
      density_ave = global_aux_var_dn%den(1)

      ! since boundary auxvar is meaningless (no pressure specified there), only use cell
#ifdef USE_ANISOTROPIC_MOBILITY
      if (dabs(dabs(dist(1))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_x
      else if (dabs(dabs(dist(2))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_y
      else if (dabs(dabs(dist(3))-1) < 1e-6) then
        ukvr = rich_aux_var_dn%kvr_z
      end if
#else
      ukvr = rich_aux_var_dn%kvr
#endif
     
      if (ukvr*perm_dn>floweps) then
        v_darcy = perm_dn * ukvr * dphi
      endif 

  end select

  q = v_darcy * area

  fluxm = q*density_ave

  Res(1)=fluxm

end subroutine RichardsBCFlux

end module Richards_Common_module
