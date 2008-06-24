module Mphase_Aux_module
use mphase_pckr_module
  implicit none
  
  private 

#include "definitions.h"
  PetscReal, parameter :: fmwnacl = 58.44277D0, Rgasj=8.31415
   

type, public :: mphase_auxvar_elem_type
   PetscReal :: pres
    PetscReal :: temp
    PetscReal , pointer :: sat(:)
    PetscReal , pointer :: den(:)
    PetscReal , pointer :: avgmw(:)
    PetscReal , pointer :: h(:)
    PetscReal , pointer :: u(:)
    PetscReal , pointer :: pc(:)
    PetscReal , pointer :: kvr(:)
    PetscReal , pointer :: xmol(:)
    PetscReal , pointer :: diff(:)
    PetscReal , pointer :: hysdat(:)
    PetscReal :: zco2
!     PetscReal :: vis
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
 end type mphase_auxvar_elem_type

  type, public :: mphase_auxvar_type
    
    type(mphase_auxvar_elem_type), pointer :: aux_var_elem(:) 
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
  end type mphase_auxvar_type
  
  type, public :: Mphase_type
     PetscInt :: n_zero_rows
     PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

     logical :: aux_vars_up_to_date
     logical :: inactive_cells_exist
     PetscInt :: num_aux, num_aux_bc
     type(Mphase_auxvar_type), pointer :: aux_vars(:)
     type(Mphase_auxvar_type), pointer :: aux_vars_bc(:)
  end type Mphase_type

  

  public :: MphaseAuxCreate, MphaseAuxDestroy, &
            MphaseAuxVarCompute_NINC, MphaseAuxVarCompute_WINC,&
            MphaseAuxVarInit, MphaseAuxVarCopy

contains
 


! ************************************************************************** !
!
! MphaseAuxVarCreate: Allocate and initialize auxilliary object
! author: Chuan Lu
! date: 02/27/08
!
! ************************************************************************** !
function MphaseAuxCreate()

  use Option_module

  implicit none
  
  type(Mphase_type), pointer :: MphaseAuxCreate
  
  type(Mphase_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = .false.
  aux%inactive_cells_exist = .false.
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  MphaseAuxCreate => aux
  
end function MphaseAuxCreate



! ************************************************************************** !
!
! MphaseAuxVarInit: Initialize auxilliary object
! author: Chuan Lu
! date: 02/14/08
!
! ************************************************************************** !
subroutine MphaseAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(Mphase_auxvar_type) :: aux_var
  type(option_type) :: option

  PetscInt :: var_elem_size, var_node_size
  PetscInt :: nvar 

  allocate(aux_var%aux_var_elem(0 : option%nflowdof))
  allocate(aux_var%aux_var_elem(0)%hysdat(4))
 
  do nvar = 0, option%nflowdof
     allocate ( aux_var%aux_var_elem(nvar)%sat(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%den(option%nphase))
     allocate ( aux_var%aux_var_elem(nvar)%avgmw(option%nphase))
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

end subroutine MphaseAuxVarInit

! ************************************************************************** !
!
! RichardsAuxVarCopy: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine MphaseAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(mphase_auxvar_elem_type) :: aux_var, aux_var2
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
  aux_var2%xmol = aux_var%xmol
  aux_var2%diff = aux_var%diff

end subroutine MphaseAuxVarCopy


! ************************************************************************** !
!
! MphaseAuxVarCompute_NI: Computes auxilliary variables for each grid cell
!                        No increments 
! author: Chuan Lu
! date: 02/22/08
!
! ************************************************************************** !
subroutine MphaseAuxVarCompute_NINC(x,aux_var,iphase,saturation_function, &
                                   fluid_properties,option)

  use Option_module
  use water_eos_module
  use gas_eos_module
  use co2eos_module
  use span_wagner_module
  use span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use Material_module
  use mphase_pckr_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(mphase_auxvar_elem_type) :: aux_var
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal ::  p, t, temp, p2, err
  PetscReal :: henry
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: xm_nacl, x_nacl, vphi             
   
  
  
  aux_var%sat = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%xmol = 0.d0
  aux_var%pc = 0.d0
  aux_var%kvr = 0.d0
  aux_var%diff = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)  
  aux_var%temp = x(2)

  p= aux_var%pres
  t= aux_var%temp
  select case(iphase)
!******* Only aqueous phase exist ***********  
    case(1)
        aux_var%xmol(2)=x(3)
        if(aux_var%xmol(2)<0.D0) print *,'tran:',iphase, x(1:3)
        if(aux_var%xmol(2)>1.D0) print *,'tran:',iphase, x(1:3)
        aux_var%xmol(1)=1.D0 - aux_var%xmol(2)
        aux_var%pc(:)=0.D0
        aux_var%sat(1)=1.D0
        aux_var%sat(2)= 0.D0
        kr(1)= 1.D0
        kr(2)= 0.D0
!******* Only gas phase exist ***********  
     case(2)
        aux_var%xmol(4)=x(3)
        if(aux_var%xmol(4)<0.D0) print *,'tran:',iphase, x(1:3)
        if(aux_var%xmol(4)>1.D0) print *,'tran:',iphase, x(1:3)
        aux_var%xmol(3)=1.D0 - aux_var%xmol(4)
        aux_var%pc(:)=0.D0
        aux_var%sat(1)= 0.D0
        aux_var%sat(2)= 1.D0
        aux_var%pc(2)=0.D0
        kr(1)= 0.D0
        kr(2)= 1.D0
    case(3)    
        aux_var%sat(2)=x(3)
        if(aux_var%sat(2)< 0.D0)then
           print *,'tran:',iphase, x(1:3)
           aux_var%sat(2)= 0.D0
        endif
        if(aux_var%sat(2)> 1.D0) print *,'tran:',iphase, x(1:3)
        aux_var%sat(1)=1.D0 - aux_var%sat(2)
        aux_var%pc(:)=0.D0
        temp = 1D-2
        aux_var%xmol(1)=1.D0; aux_var%xmol(2)=0.D0
        aux_var%xmol(3)=temp; aux_var%xmol(4)=1.D0-aux_var%xmol(3)
   end select
! ********************* Gas phase properties ***********************
    call PSAT(t, sat_pressure, ierr)
    err=1.D0
    p2 = p

    if(p2>=5d4)then
       
       if(option%co2eos == EOS_SPAN_WAGNER)then
! ************ Span-Wagner EOS ********************             
          select case(option%itable)  
          case(0,1,2,4,5)
             if( option%itable >=4) then
                ! print *,' interp', itable
                call co2_sw_interp(p2*1.D-6, t,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
             else
                call co2_span_wagner(p2*1.D-6, t +273.15D0,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
             endif
             dg= dg / option%fmwco2
             fg= fg * 1.D6 
             hg= hg * option%fmwco2
             xphi = fg/p2
! ************* Span-Wagner EOS with Bi-Cubic Spline interpolation ********
          case(3) 
             call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
             call visco2(t, dg, visg)
             dg= dg / option%fmwco2
             fg= fg * 1.D6 
             hg= hg * option%fmwco2
             xphi = fg/p2
          end select
       elseif(option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
          call CO2(t, p2,  dg,fg, xphi, hg)
          call visco2( t,dg,visg)
          dg = dg / option%fmwco2
          hg = hg * option%fmwco2 *option%scale
          !      print *, 'translator', p2, t, dg,hg,visg
       else
         print *,'pflow mphase ERROR: Need specify CO2 EOS'
         STOP    
      endif
   else      
      call ideal_gaseos_noderiv(p2, t,option%scale,dg,hg,eng)
      call visco2(t,dg*option%fmwco2,visg)
      fg=p2
      xphi = 1.D0
   endif

   call Henry_duan_sun(t, p2 *1D-5, henry, xphi, option%m_nacl, option%m_nacl,sat_pressure*1D-5)
   henry= 1D0 / (option%fmwh2o *1D-3) / (henry*1D-5 )/xphi 
   
   select case(iphase)     
   case(1)
      aux_var%xmol(4)=aux_var%xmol(2)*henry/p   
      aux_var%xmol(3)=1.D0-aux_var%xmol(4)
      if(aux_var%xmol(3)<0.D0)aux_var%xmol(3)=0.D0
     !     if(xmol(3)<0.D0) xmol(3)=0.D0
   case(2)   
      aux_var%xmol(2)= p*aux_var%xmol(4)/henry
      aux_var%xmol(1)=1.D0-aux_var%xmol(2)
   case(3)
      temp= sat_pressure / p
      aux_var%xmol(2)=(1.D0-temp)/(Henry/ p - temp)
      aux_var%xmol(1)= 1.D0- aux_var%xmol(2)
      aux_var%xmol(3)=aux_var%xmol(1) * temp
      aux_var%xmol(4)= 1.D0-aux_var%xmol(3)            
   end select
   aux_var%avgmw(2)= aux_var%xmol(3)* option%fmwh2o + aux_var%xmol(4) * option%fmwco2
   pw = p
   call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,option%scale,ierr) 
   aux_var%den(2)= 1.D0/(aux_var%xmol(4)/dg + aux_var%xmol(3)/dw_mol)
   aux_var%h(2)=  hg  
   aux_var%u(2)=  hg - p/dg * option%scale
   aux_var%pc(2)=0D0
   aux_var%diff(option%nflowspec+1:option%nflowspec*2)= 2.13D-5
!       fluid_properties%diff_base(2)
! Note: not temperature dependent yet.       
   aux_var%zco2=aux_var%den(2)/(p/rgasj/(t+273.15D0)*1D-3)
!***************  Liquid phase properties **************************
 
!  avgmw(1)= xmol(1)* fmwh2o + xmol(2) * fmwco2 
  aux_var%h(1) = hw
  aux_var%u(1) = aux_var%h(1) - pw /dw_mol* option%scale
  aux_var%diff(1:option%nflowspec) = 1D-9
  ! fluid_properties%diff_base(1)

  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.D3 + xm_nacl)
  call nacl_den(t, p*1D-6, xm_nacl, dw_kg) 
  dw_kg = dw_kg * 1D3
  call nacl_vis(t,p*1D-6,xm_nacl,visl)

!FEHM mixing ****************************
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
  !den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing **************************
  x_nacl =  option%m_nacl/( option%m_nacl + 1D3/option%fmwh2o)
! **  xmol(1) = xh2o + xnacl
  aux_var%avgmw(1)= (aux_var%xmol(1) - x_nacl) * option%fmwh2o&
       + x_nacl * fmwnacl + aux_var%xmol(2) * option%fmwco2
  vphi=1D-6*(37.51D0 + t&
       *(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))
  aux_var%den(1)=dw_kg/(1D0-(option%fmwco2*1D-3-dw_kg*vphi)&
       *aux_var%xmol(2)/(aux_var%avgmw(1)*1D-3))
  aux_var%den(1)=aux_var%den(1)/aux_var%avgmw(1)
  
       
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 ***********
 !   den(1)= 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1)=dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1)=den(1)/avgmw(1)
!******************************** 2 phase S-Pc-kr relation ***********************************
    if(option%nphase>=2)then
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
   endif

!    call SaturationFunctionCompute(aux_var%pres,aux_var%sat,kr, &
!                                   ds_dp,dkr_dp, &
!                                   saturation_function, &
!                                   por,perm, &
!                                   option)
       aux_var%kvr(2)=kr(2)/visg     
       aux_var%kvr(1) = kr(1)/visl
       select case(iphase)
         case(1)
           aux_var%pc =0.D0
         case(2)
           aux_var%pc =0.D0
      end select          

    end subroutine MphaseAuxVarCompute_NINC



subroutine MphaseAuxVarCompute_WINC(x, delx, aux_var,iphase,saturation_function, &
                                    fluid_properties,option)

  use Option_module
  use water_eos_module
  use Material_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  type(mphase_auxvar_elem_type) :: aux_var(1:option%nflowdof)
  PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  MphaseAuxVarCompute_NINC(xx,aux_var(n),iphase,saturation_function, &
                                   fluid_properties, option)
  enddo

end subroutine MphaseAuxVarCompute_WINC

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MphaseAuxVarDestroy(aux_var)

  implicit none

  type(mphase_auxvar_elem_type) :: aux_var
  
  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
  nullify(aux_var%xmol)
  if (associated(aux_var%diff))deallocate(aux_var%diff)
  nullify(aux_var%diff)
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
 if (associated(aux_var%avgmw))deallocate(aux_var%avgmw)
  nullify(aux_var%u)
end subroutine MphaseAuxVarDestroy

! ************************************************************************** !
!
! RichardsAuxDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MphaseAuxDestroy(aux, option)
use option_module
  implicit none

  type(mphase_type), pointer :: aux
  type(option_type) :: option
  PetscInt :: iaux, ielem
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    do ielem= 0, option%nflowdof 
      call MphaseAuxVarDestroy(aux%aux_vars(iaux)%aux_var_elem(ielem))
    enddo
  enddo  
  do iaux = 1, aux%num_aux_bc
    do ielem= 0, option%nflowdof 
      call MphaseAuxVarDestroy(aux%aux_vars_bc(iaux)%aux_var_elem(ielem))
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
    
end subroutine MphaseAuxDestroy



end module Mphase_Aux_module


