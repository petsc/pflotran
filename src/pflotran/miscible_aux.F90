module Miscible_Aux_module
  implicit none
  
  private 

#include "definitions.h"

type, public :: Miscible_auxvar_elem_type
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
 end type Miscible_auxvar_elem_type

  type, public :: Miscible_auxvar_type
    
    type(Miscible_auxvar_elem_type), pointer :: aux_var_elem(:) 
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
  end type Miscible_auxvar_type
  
  type, public :: Miscible_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: ckdry(:)
    PetscReal, pointer :: sir(:,:)
  end type Miscible_parameter_type
    
  type, public :: Miscible_type
     PetscInt :: n_zero_rows
     PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

     PetscBool :: aux_vars_up_to_date
     PetscBool :: inactive_cells_exist
     PetscInt :: num_aux, num_aux_bc
     type(Miscible_parameter_type), pointer :: Miscible_parameter
     type(Miscible_auxvar_type), pointer :: aux_vars(:)
     type(Miscible_auxvar_type), pointer :: aux_vars_bc(:)
     PetscReal , pointer :: Resold_AR(:,:)
     PetscReal , pointer :: Resold_BC(:,:)
     PetscReal , pointer :: Resold_FL(:,:)
     PetscReal , pointer :: delx(:,:)
  end type Miscible_type

  

  public :: MiscibleAuxCreate, MiscibleAuxDestroy, &
            MiscibleAuxVarCompute_NINC, MiscibleAuxVarCompute_WINC,&
            MiscibleAuxVarInit, MiscibleAuxVarCopy

contains
 


! ************************************************************************** !
!
! MiscibleAuxVarCreate: Allocate and initialize auxilliary object
! author: Chuan Lu
! date: 02/27/08
!
! ************************************************************************** !
function MiscibleAuxCreate()

  use Option_module

  implicit none
  
  type(Miscible_type), pointer :: MiscibleAuxCreate
  
  type(Miscible_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  allocate(aux%Miscible_parameter)
  nullify(aux%Miscible_parameter%sir)
  nullify(aux%Miscible_parameter%ckwet)
  nullify(aux%Miscible_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  MiscibleAuxCreate => aux
  
end function MiscibleAuxCreate



! ************************************************************************** !
!
! MiscibleAuxVarInit: Initialize auxilliary object
! author: Chuan Lu
! date: 02/14/08
!
! ************************************************************************** !
subroutine MiscibleAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(Miscible_auxvar_type) :: aux_var
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

end subroutine MiscibleAuxVarInit

! ************************************************************************** !
!
! MiscibleAuxVarCopy: Copies an auxilliary variable
! author: Chuan Lu
! date: 10/13/0
!
! ************************************************************************** !  
subroutine MiscibleAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(Miscible_auxvar_elem_type) :: aux_var, aux_var2
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

end subroutine MiscibleAuxVarCopy


subroutine Water_glycal_density( y,p, dkg)
  implicit none
  PetscReal y, p ! water mass fraction
  PetscReal dkg

  dkg = ((0.0806*y-0.203)*y + 0.0873)*y + 1.0341D0
  dkg = (4.49758D-10* y +(1D0-y)*5D-10)*(p-1.01325D5) + dkg
  dkg = dkg * 1D3  ! convert g/cm^3 to kg/m^3
 
end subroutine Water_glycal_density       
! ************************************************************************** !
!
! MiscibleAuxVarCompute_NI: Computes auxilliary variables for each grid cell
!                        No increments 
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleAuxVarCompute_NINC(x,aux_var,global_aux_var, &
             fluid_properties,option)

  use Option_module
  use Global_Aux_module  
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
!  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(Miscible_auxvar_elem_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscInt :: iphase
 

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal ::  p, t, temp, p2, err
  PetscReal :: henry,lngamco2
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: tk, visw 
  PetscReal :: denw, yh2o
  PetscReal :: tmp 
  PetscInt :: iflag  
  
  aux_var%sat = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%pc = 0.d0
  aux_var%kvr = 0.d0
  aux_var%xmol = 0.d0
! aux_var%diff = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)  
  aux_var%xmol(2:option%nflowspec) = x(2:option%nflowspec)
  tmp=sum(aux_var%xmol)
  aux_var%xmol(1) = 1D0 - tmp

  aux_var%avgmw(1) = aux_var%xmol(1)*FMWH2O + aux_var%xmol(2)*FMWGLYC
  yh2o = aux_var%xmol(1)*FMWH2O/aux_var%avgmw(1)
  call Water_glycal_density(yh2o, aux_var%pres, denw)
  aux_var%den(1)=denw/aux_var%avgmw(1)
  visw = 10D0**(1.6743d0*yh2o-0.0758) * 1.0D-3
  aux_var%sat(1)=1D0
  aux_var%kvr(1) = 1.d0/visw
  aux_var%h(1)= denw * 4.18D-3*global_aux_var%temp(1)
!  auc_var%fdiff(1) = ((((-4.021*y+9.1181)*y-5.9703)*y+0.4043D-3)*y + 0.5687)*1D-5
  aux_var%diff(1:option%nflowspec) =  fluid_properties%diffusion_coefficient
  

end subroutine MiscibleAuxVarCompute_NINC



subroutine MiscibleAuxVarCompute_WINC(x, delx, aux_var,global_auxvar,&
                                    fluid_properties,option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
!  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  type(Miscible_auxvar_elem_type) :: aux_var(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar
 ! PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  MiscibleAuxVarCompute_NINC(xx,aux_var(n),global_auxvar, &
      fluid_properties, option)
  enddo

end subroutine MiscibleAuxVarCompute_WINC

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a Miscible auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MiscibleAuxVarDestroy(aux_var)

  implicit none

  type(Miscible_auxvar_elem_type) :: aux_var
  
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
end subroutine MiscibleAuxVarDestroy

! ************************************************************************** !
!
! RichardsAuxDestroy: Deallocates a Miscible auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine MiscibleAuxDestroy(aux, option)
use option_module
  implicit none

  type(Miscible_type), pointer :: aux
  type(option_type) :: option
  PetscInt :: iaux, ielem
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    do ielem= 0, option%nflowdof 
      call MiscibleAuxVarDestroy(aux%aux_vars(iaux)%aux_var_elem(ielem))
    enddo
  enddo  
  do iaux = 1, aux%num_aux_bc
    do ielem= 0, option%nflowdof 
      call MiscibleAuxVarDestroy(aux%aux_vars_bc(iaux)%aux_var_elem(ielem))
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
  if (associated(aux%Miscible_parameter)) then
    if (associated(aux%Miscible_parameter%dencpr)) deallocate(aux%Miscible_parameter%dencpr)
    nullify(aux%Miscible_parameter%dencpr)
    if (associated(aux%Miscible_parameter%ckwet)) deallocate(aux%Miscible_parameter%ckwet)
    nullify(aux%Miscible_parameter%ckwet)
    if (associated(aux%Miscible_parameter%ckdry)) deallocate(aux%Miscible_parameter%ckdry)
    nullify(aux%Miscible_parameter%ckdry)
    if (associated(aux%Miscible_parameter%sir)) deallocate(aux%Miscible_parameter%sir)
    nullify(aux%Miscible_parameter%sir)
    deallocate(aux%Miscible_parameter)
  endif
  nullify(aux%Miscible_parameter%dencpr)
    
end subroutine MiscibleAuxDestroy



end module Miscible_Aux_module


