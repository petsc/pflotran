module Flash2_Aux_module
use Mphase_pckr_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 
!#define GARCIA 1
#define DUANDEN 1

#include "finclude/petscsys.h"

type, public :: Flash2_auxvar_elem_type
   PetscReal :: pres
    PetscReal :: temp
    PetscReal , pointer :: sat(:)
    PetscReal , pointer :: den(:)
    PetscReal , pointer :: avgmw(:)
    PetscReal , pointer :: vis(:)
    PetscReal , pointer :: h(:)
    PetscReal , pointer :: u(:)
    PetscReal , pointer :: pc(:)
    PetscReal , pointer :: kvr(:)
    PetscReal , pointer :: xmol(:)
    PetscReal , pointer :: diff(:)
    PetscReal , pointer :: hysdat(:)
    PetscReal :: zco2
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
 end type Flash2_auxvar_elem_type

  type, public :: Flash2_auxvar_type
    
    type(Flash2_auxvar_elem_type), pointer :: aux_var_elem(:) 
#if 0
    PetscReal , pointer :: davgmw_dx(:)
    PetscReal , pointer :: dden_dp(:)
    PetscReal , pointer :: dden_dt(:)
    PetscReal , pointer :: dden_dx(:)
    PetscReal , pointer :: dkvr_dp(:)
    PetscReal , pointer :: dkvr_dt(:)
    PetscReal , pointer :: dkvr_ds(:)
    PetscReal , pointer :: dkvr_dx(:)
    PetscReal , pointer :: dh_dp(:)
    PetscReal , pointer :: dh_dt(:)
    PetscReal , pointer :: dh_dx(:)
    PetscReal , pointer :: du_dp(:)
    PetscReal , pointer :: du_dt(:)
    PetscReal , pointer :: du_dx(:)
#endif
  end type Flash2_auxvar_type
  
  type, public :: Flash2_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: ckdry(:)
    PetscReal, pointer :: sir(:,:)
  end type Flash2_parameter_type
    
  type, public :: Flash2_type
     PetscInt :: n_zero_rows
     PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

     PetscBool :: aux_vars_up_to_date
     PetscBool :: inactive_cells_exist
     PetscInt :: num_aux, num_aux_bc
     type(Flash2_parameter_type), pointer :: Flash2_parameter
     type(Flash2_auxvar_type), pointer :: aux_vars(:)
     type(Flash2_auxvar_type), pointer :: aux_vars_bc(:)
     PetscReal , pointer :: Resold_AR(:,:)
     PetscReal , pointer :: Resold_BC(:,:)
     PetscReal , pointer :: Resold_FL(:,:)
     PetscReal , pointer :: delx(:,:)
  end type Flash2_type

  

  public :: Flash2AuxCreate, Flash2AuxDestroy, &
            Flash2AuxVarCompute_NINC, Flash2AuxVarCompute_WINC,&
            Flash2AuxVarInit, Flash2AuxVarCopy

contains
 


! ************************************************************************** !
!
! Flash2AuxVarCreate: Allocate and initialize auxiliary object
! author: Chuan Lu
! date: 02/27/08
!
! ************************************************************************** !
function Flash2AuxCreate()

  use Option_module

  implicit none
  
  type(Flash2_type), pointer :: Flash2AuxCreate
  
  type(Flash2_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  allocate(aux%Flash2_parameter)
  nullify(aux%Flash2_parameter%sir)
  nullify(aux%Flash2_parameter%ckwet)
  nullify(aux%Flash2_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  Flash2AuxCreate => aux
  
end function Flash2AuxCreate



! ************************************************************************** !
!
! Flash2AuxVarInit: Initialize auxiliary object
! author: Chuan Lu
! date: 02/14/08
!
! ************************************************************************** !
subroutine Flash2AuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(Flash2_auxvar_type) :: aux_var
  type(option_type) :: option

  PetscInt :: var_elem_size, var_node_size
  PetscInt :: nvar 

  allocate(aux_var%aux_var_elem(0 : option%nflowdof))
  allocate(aux_var%aux_var_elem(0)%hysdat(4))
 
  do nvar = 0, option%nflowdof
     allocate ( aux_var%aux_var_elem(nvar)%sat(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%den(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%avgmw(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%vis(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%h(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%u(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%pc(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%kvr(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%xmol(option%nphase*option%nflowspec))
     allocate ( aux_var%aux_var_elem(nvar)%diff(option%nphase*option%nflowspec))
     if(nvar>0)&
     aux_var%aux_var_elem(nvar)%hysdat => aux_var%aux_var_elem(0)%hysdat

     aux_var%aux_var_elem(nvar)%pres = 0.d0
     aux_var%aux_var_elem(nvar)%temp = 0.d0
     aux_var%aux_var_elem(nvar)%sat = 0.d0
     aux_var%aux_var_elem(nvar)%den = 0.d0
     aux_var%aux_var_elem(nvar)%avgmw = 0.d0
     aux_var%aux_var_elem(nvar)%vis = 0.d0
     aux_var%aux_var_elem(nvar)%h = 0.d0
     aux_var%aux_var_elem(nvar)%u = 0.d0
     aux_var%aux_var_elem(nvar)%pc = 0.d0
     aux_var%aux_var_elem(nvar)%kvr = 0.d0
     aux_var%aux_var_elem(nvar)%xmol = 0.d0
     aux_var%aux_var_elem(nvar)%diff = 0.d0
#if 0
     aux_var%aux_var_elem(nvar)%dsat_dp = 0.d0
     aux_var%aux_var_elem(nvar)%dden_dp = 0.d0
     aux_var%aux_var_elem(nvar)%dkvr_dp = 0.d0
#endif
  enddo

end subroutine Flash2AuxVarInit

! ************************************************************************** !
!
! Flash2AuxVarCopy: Copies an auxiliary variable
! author: Chuan Lu
! date: 10/13/0
!
! ************************************************************************** !  
subroutine Flash2AuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(Flash2_auxvar_elem_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
!  aux_var2%kr = aux_var%kr
!  aux_var2%dkr_dp = aux_var%dkr_dp
!  aux_var2%vis = aux_var%vis
!  aux_var2%dvis_dp = aux_var%dvis_dp
  aux_var2%kvr = aux_var%kvr
#if 0
  aux_var2%dsat_dp = aux_var%dsat_dp
  aux_var2%dden_dp = aux_var%dden_dp
  aux_var2%dden_dt = aux_var%dden_dt
  aux_var2%dkvr_dp = aux_var%dkvr_dp
  aux_var2%dkvr_dt = aux_var%dkvr_dt
  aux_var2%dh_dp = aux_var%dh_dp
  aux_var2%dh_dt = aux_var%dh_dt
  aux_var2%du_dp = aux_var%du_dp
  aux_var2%du_dt = aux_var%du_dt  
#endif
!  aux_var2%xmol = aux_var%xmol
!  aux_var2%diff = aux_var%diff

end subroutine Flash2AuxVarCopy


! ************************************************************************** !
!
! Flash2AuxVarCompute_NI: Computes auxiliary variables for each grid cell
!                        No increments 
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine Flash2AuxVarCompute_NINC(x,aux_var,global_aux_var, &
             saturation_function,fluid_properties,option,xphico2)

  use Option_module
  use Global_Aux_module  
  use EOS_Water_module
  use Gas_EOS_module
  use co2eos_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use Saturation_Function_module
  use Fluid_module
  use Mphase_pckr_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(Flash2_auxvar_elem_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscInt :: iphase
  PetscReal, optional :: xphico2

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal ::  p, t, temp, p2, err
  PetscReal :: henry,lngamco2
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: m_na,m_cl,m_nacl, xm_nacl, x_nacl, y_nacl, vphi             
  PetscReal :: tk, xco2, pw_kg, x1, vphi_a1, vphi_a2 
  PetscReal :: Qkco2, mco2,xco2eq
  PetscReal :: tmp 
  PetscInt :: iflag  
  
  aux_var%sat = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%pc = 0.d0
  aux_var%kvr = 0.d0
! aux_var%xmol = 0.d0
! aux_var%diff = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)  
  aux_var%temp = x(2)

  p = aux_var%pres
  t = aux_var%temp

! ********************* Gas phase properties ***********************
    call EOSWaterSaturationPressure(t, sat_pressure, ierr)
    err=1.D0
    p2 = p

    if(p2 >= 5.d4)then
      if(option%co2eos == EOS_SPAN_WAGNER)then
! ************ Span-Wagner EOS ********************             
        select case(option%itable)  
          case(0,1,2,4,5)
            if( option%itable >=4) then
                ! print *,' interp', itable
              call co2_sw_interp(p2*1.D-6, t,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            else
              iflag = 1
              call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag, &
                     option%itable)
            endif
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
            xphi = fg/p2
! ************* Span-Wagner EOS with Bi-Cubic Spline interpolation ********
          case(3) 
            call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
            call visco2(t, dg, visg)
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
            xphi = fg/p2
        end select
      elseif(option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
        call CO2(t, p2,  dg,fg, xphi, hg)
        call visco2( t,dg,visg)
        dg = dg / FMWCO2
        hg = hg * FMWCO2 *option%scale
          !      print *, 'translator', p2, t, dg,hg,visg
      else
        call printErrMsg(option,'pflow Flash2 ERROR: Need specify CO2 EOS')
      endif
    else      
      call ideal_gaseos_noderiv(p2, t,option%scale,dg,hg,eng)
      call visco2(t,dg*FMWCO2,visg)
      fg=p2
      xphi = 1.D0
    endif
 
!*********** Get Salinity properties ***********************
    m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na 
    if (option%ntrandof>0) then
      m_na = global_aux_var%m_nacl(1)
      m_cl = global_aux_var%m_nacl(2)
      m_nacl = m_na
      if (m_cl> m_na) m_nacl = m_cl
    endif  

!** Calculate solubility of CO2 in aqueoues phase *************
    call Henry_duan_sun(t, p2 *1D-5, henry, xphi, lngamco2, &
           m_na,m_cl,sat_pressure*1D-5)

    Qkco2 = henry*xphi  ! convert from bar to Pa
    henry = 1.D0 / (FMWH2O*1.D-3) / (henry*1.D-5) / xphi 
    if(present(xphico2)) xphico2 = xphi
   
    mco2 = (p - sat_pressure)*1D-5 * Qkco2
    xco2eq = mco2/(1D3/fmwh2o + mco2 + m_nacl) 
   
    tmp= Henry/p
    if (x(3) < xco2eq) then
      ! water only
      aux_var%xmol(2)=x(3)
      aux_var%xmol(1)=1.D0 - aux_var%xmol(2)
      aux_var%xmol(4)=aux_var%xmol(2)*tmp
      aux_var%xmol(3)=1.D0 - aux_var%xmol(4)  
      aux_var%sat(1)=1.D0
      aux_var%sat(2)=0.D0
      iphase = 1
    elseif (x(3) > (1.D0-sat_pressure/p)) then
	    !gas only
      iphase =2
      aux_var%xmol(4)=x(3)
      aux_var%xmol(3)=1.D0 - aux_var%xmol(4) 
      aux_var%xmol(2)=aux_var%xmol(4)/tmp
      aux_var%xmol(1)=1.D0 - aux_var%xmol(2)
      aux_var%sat(1)=0.D0 !1.D-8
      aux_var%sat(2)=1.D0
    else 
      iphase = 3
      aux_var%xmol(1)=1.D0 - xco2eq
      aux_var%xmol(2)= xco2eq
      aux_var%xmol(3)= sat_pressure/p*aux_var%xmol(1) 
      aux_var%xmol(4)= 1.D0 - aux_var%xmol(3)
    endif 

! **************  Gas phase properties ********************
    aux_var%avgmw(2) = aux_var%xmol(3)*FMWH2O + aux_var%xmol(4)*FMWCO2
    pw = p
    call EOSWaterDensityEnthalpy(t,pw,dw_kg,dw_mol,hw,option%scale,ierr)
    aux_var%den(2) = 1.D0/(aux_var%xmol(4)/dg + aux_var%xmol(3)/dw_mol)
    aux_var%h(2) = hg  
    aux_var%u(2) = hg - p/dg * option%scale
    
!   aux_var%diff(option%nflowspec+1:option%nflowspec*2) = 2.13D-5
    aux_var%diff(option%nflowspec+1:option%nflowspec*2) = &
      fluid_properties%gas_diffusion_coefficient
      
!       fluid_properties%diff_base(2)
! Note: not temperature dependent yet.       
!  z factor    
    aux_var%zco2=aux_var%den(2)/(p/IDEAL_GAS_CONST/(t+273.15D0)*1D-3)

 !***************  Liquid phase properties **************************
 
!    avgmw(1)= xmol(1)* FMWH2O + xmol(2) * FMWCO2 
    aux_var%h(1) = hw
    aux_var%u(1) = aux_var%h(1) - pw /dw_mol* option%scale
    
    aux_var%diff(1:option%nflowspec) = fluid_properties%diffusion_coefficient
  ! fluid_properties%diff_base(1) need more work here. Add temp. dependence.
  
    xm_nacl = m_nacl * FMWNACL
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)
    call EOSWaterDensityNaCl(t, p*1D-6, xm_nacl, dw_kg) 
    dw_kg = dw_kg * 1D3
    call EOSWaterViscosityNaCl(t,p*1D-6,xm_nacl,visl)
    
    y_nacl =  m_nacl/( m_nacl + 1D3/FMWH2O)
!   y_nacl is the mole fraction
    aux_var%avgmw(1)= aux_var%xmol(1)*((1D0 - y_nacl) * FMWH2O&
       + y_nacl * FMWNACL) + aux_var%xmol(2) * FMWCO2

!duan mixing **************************
#ifdef DUANDEN
  call EOSWaterDuanMixture (t,p,aux_var%xmol(2),y_nacl,aux_var%avgmw(1),dw_kg,aux_var%den(1))
#endif 

! Garcia mixing **************************
#ifdef GARCIA
  vphi=1D-6*(37.51D0 + t&
       *(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))
  aux_var%den(1)=dw_kg/(1D0-(FMWCO2*1D-3-dw_kg*vphi)&
       *aux_var%xmol(2)/(aux_var%avgmw(1)*1D-3))
  aux_var%den(1)=aux_var%den(1)/aux_var%avgmw(1)
#endif  
       

!FEHM mixing ****************************
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
  !den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing **************************
        
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 ***********
 !   den(1)= 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1)=dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1)=den(1)/avgmw(1)

!****** calcultate phase splition for 2 phase coexist condition *******
  select case(iphase)
  case(1,2)
  case(3)
    aux_var%sat(2) = aux_var%den(1)* ( x(3) - aux_var%xmol(2))/&
      (aux_var%den(2) * (aux_var%xmol(4)-x(3)) - aux_var%den(1)*(aux_var%xmol(2)-x(3)))
    if(aux_var%sat(2) >1D0 .or. aux_var%sat(2) <0D0) print *,'z->s error: ',aux_var%sat(2)
    if(aux_var%sat(2) > 1D0) aux_var%sat(2) = 1D0
    if(aux_var%sat(2) < 0D0) aux_var%sat(2) = 0D0  
    aux_var%sat(1) = 1D0 - aux_var%sat(2)
  end select
 
!******************************** 2 phase S-Pc-kr relation ********************
    aux_var%pc =0.D0

      if(saturation_function%hysteresis_id <=0.1D0 ) then 
         call pckrNH_noderiv(aux_var%sat,aux_var%pc,kr, &
                                   saturation_function, &
                                   option)
        pw=p !-pc(1)
     
       else
          call pckrHY_noderiv(aux_var%sat,aux_var%hysdat,aux_var%pc,kr, &
                                   saturation_function, &
                                   option)
     end if

!    call SaturationFunctionCompute(aux_var%pres,aux_var%sat,kr, &
!                                   ds_dp,dkr_dp, &
!                                   saturation_function, &
!                                   por,perm, &
!                                   option)
       aux_var%kvr(2) = kr(2)/visg     
       aux_var%kvr(1) = kr(1)/visl
       aux_var%vis(2) = visg     
       aux_var%vis(1) = visl

end subroutine Flash2AuxVarCompute_NINC



subroutine Flash2AuxVarCompute_WINC(x, delx, aux_var,global_auxvar,saturation_function, &
                                    fluid_properties,option)

  use Option_module
  use Global_Aux_module
  
  use Saturation_Function_module
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  type(Flash2_auxvar_elem_type) :: aux_var(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar
 ! PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  Flash2AuxVarCompute_NINC(xx,aux_var(n),global_auxvar, &
      saturation_function,fluid_properties, option)
  enddo

end subroutine Flash2AuxVarCompute_WINC

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a FLASH2 auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine Flash2AuxVarDestroy(aux_var)

  implicit none

  type(Flash2_auxvar_elem_type) :: aux_var
  
!  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
!  nullify(aux_var%xmol)
!  if (associated(aux_var%diff))deallocate(aux_var%diff)
!  nullify(aux_var%diff)
  if (associated(aux_var%pc))deallocate(aux_var%pc)
  nullify(aux_var%pc)
  if (associated(aux_var%sat))deallocate(aux_var%sat)
  nullify(aux_var%sat)
  if (associated(aux_var%u))deallocate(aux_var%u)
  nullify(aux_var%u)
  if (associated(aux_var%h))deallocate(aux_var%h)
  nullify(aux_var%h)
  if (associated(aux_var%den))deallocate(aux_var%den)
  nullify(aux_var%den)
  if (associated(aux_var%den))deallocate(aux_var%vis)
  nullify(aux_var%vis)
  if (associated(aux_var%avgmw))deallocate(aux_var%avgmw)
  nullify(aux_var%avgmw)
end subroutine Flash2AuxVarDestroy

! ************************************************************************** !
!
! RichardsAuxDestroy: Deallocates a FLASH2 auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine Flash2AuxDestroy(aux, option)

  use Option_module
  implicit none

  type(Flash2_type), pointer :: aux
  type(option_type) :: option
  PetscInt :: iaux, ielem
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    do ielem= 0, option%nflowdof 
      call Flash2AuxVarDestroy(aux%aux_vars(iaux)%aux_var_elem(ielem))
    enddo
  enddo  
  do iaux = 1, aux%num_aux_bc
    do ielem= 0, option%nflowdof 
      call Flash2AuxVarDestroy(aux%aux_vars_bc(iaux)%aux_var_elem(ielem))
    enddo
  enddo  
  
  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) deallocate(aux%aux_vars_bc)
  nullify(aux%aux_vars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%Flash2_parameter)) then
    if (associated(aux%Flash2_parameter%dencpr)) deallocate(aux%Flash2_parameter%dencpr)
    nullify(aux%Flash2_parameter%dencpr)
    if (associated(aux%Flash2_parameter%ckwet)) deallocate(aux%Flash2_parameter%ckwet)
    nullify(aux%Flash2_parameter%ckwet)
    if (associated(aux%Flash2_parameter%ckdry)) deallocate(aux%Flash2_parameter%ckdry)
    nullify(aux%Flash2_parameter%ckdry)
    if (associated(aux%Flash2_parameter%sir)) deallocate(aux%Flash2_parameter%sir)
    nullify(aux%Flash2_parameter%sir)
    deallocate(aux%Flash2_parameter)
  endif
  nullify(aux%Flash2_parameter%dencpr)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine Flash2AuxDestroy



end module Flash2_Aux_module


