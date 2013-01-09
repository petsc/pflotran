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
 end type Miscible_auxvar_elem_type

  type, public :: Miscible_auxvar_type
    type(Miscible_auxvar_elem_type), pointer :: aux_var_elem(:) 
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
     PetscInt :: num_aux, num_aux_bc, num_aux_ss
     type(Miscible_parameter_type), pointer :: Miscible_parameter
     type(Miscible_auxvar_type), pointer :: aux_vars(:)
     type(Miscible_auxvar_type), pointer :: aux_vars_bc(:)
     type(Miscible_auxvar_type), pointer :: aux_vars_ss(:)
     PetscReal, pointer :: Resold_AR(:,:)
     PetscReal, pointer :: Resold_BC(:,:)
     PetscReal, pointer :: Resold_FL(:,:)
     PetscReal, pointer :: delx(:,:)
  end type Miscible_type

  public :: MiscibleAuxCreate, MiscibleAuxDestroy, &
            MiscibleAuxVarCompute_NINC, MiscibleAuxVarCompute_WINC,&
            MiscibleAuxVarInit, MiscibleAuxVarCopy

contains

! ************************************************************************** !
!
! MiscibleAuxVarCreate: Allocate and initialize auxiliary object
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
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
  aux%n_zero_rows = 0
  allocate(aux%Miscible_parameter)
  nullify(aux%Miscible_parameter%sir)
  nullify(aux%Miscible_parameter%ckwet)
  nullify(aux%Miscible_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
  nullify(aux%Resold_AR)
  nullify(aux%Resold_BC)
  nullify(aux%Resold_FL)
  nullify(aux%delx)

  MiscibleAuxCreate => aux
  
end function MiscibleAuxCreate



! ************************************************************************** !
!
! MiscibleAuxVarInit: Initialize auxiliary object
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
  enddo

end subroutine MiscibleAuxVarInit

! ************************************************************************** !
!
! MiscibleAuxVarCopy: Copies an auxiliary variable
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
  aux_var2%kvr = aux_var%kvr
!  aux_var2%xmol = aux_var%xmol
!  aux_var2%diff = aux_var%diff

end subroutine MiscibleAuxVarCopy

! ************************************************************************** !
!
! Water_glycol_density: Computes water-propylene glycol mixture density 
! author: Chuan Lu
! date: 12/12/11
!
! ************************************************************************** !
subroutine Water_glycol_density(y,p,dkg)
  implicit none
  PetscReal y, p ! water mass fraction
  PetscReal dkg

  dkg = (((0.0806d0*y-0.203d0)*y + 0.0873d0)*y + 1.0341d0) * 1.d3
! dkg = (4.49758d-10*y +(1.d0-y)*5.d-10)*(p-1.01325d5) + dkg
! dkg = dkg * 1.d3  ! convert g/cm^3 to kg/m^3

  dkg = dkg * (1+(4.49758d-10*y +(1.d0-y)*5.d-10)*(p-1.01325d5))
 
end subroutine Water_glycol_density

! ************************************************************************** !
!
! MiscibleAuxVarCompute_NINC: Computes auxiliary variables for each grid cell
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
  type(Miscible_auxvar_elem_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  
  PetscReal :: x(option%nflowdof)
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: p, t, temp, p2, err
  PetscReal :: henry,lngamco2
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: tk, visw 
  PetscReal :: denw, yh2o, yppg
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

! aux_var%xmol(2) = x(2)
  aux_var%xmol(2) = exp(x(2))
! aux_var%xmol(2) = (atan(x(2))/3.1416*2+1)/2
  aux_var%xmol(1) = 1.D0 - aux_var%xmol(2)

! Glycol-Water mixture formula weight (kg/kmol)
  aux_var%avgmw(1) = aux_var%xmol(1)*FMWH2O + aux_var%xmol(2)*FMWGLYC
  
! Mass fraction water
  yh2o = aux_var%xmol(1)*FMWH2O/aux_var%avgmw(1)
  
  call Water_glycol_density(yh2o, aux_var%pres, denw)
  
  aux_var%den(1) = denw/aux_var%avgmw(1)
  
! Glycol-Water mixture viscosity (yh2o mass fraction water)
  yppg = 1.d0-yh2o
  visw = 10.d0**(1.6743d0*yppg-0.0758d0) * 1.0d-3 ! centipoise to Pa s.

  aux_var%vis(1) = visw
  
  aux_var%sat(1) = 1.d0
  aux_var%kvr(1) = 1.d0/visw
  aux_var%h(1) = denw*4.18d-3*global_aux_var%temp(1)
  
! Glycol-Water mixture diffusivity (yh2o mass fraction water)
  aux_var%diff(2) = ((((-4.021d0*yh2o + 9.1181d0)*yh2o - 5.9703d0)*yh2o &
     + 0.4043d0)*yh2o + 0.5687d0) * 1.d-9 ! m^2/s
  aux_var%diff(1) = aux_var%diff(2)

! aux_var%diff(1:option%nflowspec) = fluid_properties%diffusion_coefficient

end subroutine MiscibleAuxVarCompute_NINC



subroutine MiscibleAuxVarCompute_WINC(x,delx,aux_var,global_auxvar, &
                                    fluid_properties,option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(Miscible_auxvar_elem_type) :: aux_var(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar

  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  PetscInt :: idof 
  
  do idof = 1, option%nflowdof
    xx = x; xx(idof) = x(idof) + delx(idof)

!   print *,'Winc: ',idof,x(idof),xx(idof),delx(idof)

! ***   note: var_node here starts from 1 to option%flowdof ***
    call  MiscibleAuxVarCompute_NINC(xx,aux_var(idof),global_auxvar, &
      fluid_properties,option)
  enddo

end subroutine MiscibleAuxVarCompute_WINC


! ************************************************************************** !
!
! MiscibleAuxVarElemDestroy: Deallocates a mphase auxiliary elment object
! author: 
! date: 
!
! ************************************************************************** !
subroutine MiscibleAuxVarElemDestroy(aux_var_elem)

  implicit none

  type(miscible_auxvar_elem_type) :: aux_var_elem

  if (associated(aux_var_elem%xmol)) deallocate(aux_var_elem%xmol)
  nullify(aux_var_elem%xmol)
  if (associated(aux_var_elem%diff)) deallocate(aux_var_elem%diff)
  nullify(aux_var_elem%diff)
  if (associated(aux_var_elem%pc)) deallocate(aux_var_elem%pc)
  nullify(aux_var_elem%pc)
  if (associated(aux_var_elem%sat)) deallocate(aux_var_elem%sat)
  nullify(aux_var_elem%sat)
  if (associated(aux_var_elem%u)) deallocate(aux_var_elem%u)
  nullify(aux_var_elem%u)
  if (associated(aux_var_elem%h)) deallocate(aux_var_elem%h)
  nullify(aux_var_elem%h)
  if (associated(aux_var_elem%den)) deallocate(aux_var_elem%den)
  nullify(aux_var_elem%den)
  if (associated(aux_var_elem%den)) deallocate(aux_var_elem%vis)
  nullify(aux_var_elem%vis)
  if (associated(aux_var_elem%avgmw)) deallocate(aux_var_elem%avgmw)
  nullify(aux_var_elem%avgmw)

end subroutine MiscibleAuxVarElemDestroy

! ************************************************************************** !
!
! MiscibleAuxVarDestroy: Deallocates a miscible auxiliary object
! author: 
! date: 
!
! ************************************************************************** !
subroutine MiscibleAuxVarDestroy(aux_var)

  implicit none

  type(miscible_auxvar_type) :: aux_var

  PetscInt :: ielem

  ! subtract 1 since indexing from 0
  if (associated(aux_var%aux_var_elem)) then
    do ielem = 0, size(aux_var%aux_var_elem) - 1 
      call MiscibleAuxVarElemDestroy(aux_var%aux_var_elem(ielem))
    enddo
    deallocate(aux_var%aux_var_elem)
    nullify(aux_var%aux_var_elem)
  endif

end subroutine MiscibleAuxVarDestroy

! ************************************************************************** !
!
! MiscibleAuxDestroy: Deallocates a miscible auxiliary object
! author: 
! date: 
!
! ************************************************************************** !
subroutine MiscibleAuxDestroy(aux)

  implicit none

  type(miscible_type), pointer :: aux

  PetscInt :: iaux
  
  if (.not.associated(aux)) return

  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call MiscibleAuxVarDestroy(aux%aux_vars(iaux))
    enddo
    deallocate(aux%aux_vars)
    nullify(aux%aux_vars)
  endif
  
  if (associated(aux%aux_vars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call MiscibleAuxVarDestroy(aux%aux_vars_bc(iaux))
    enddo
    deallocate(aux%aux_vars_bc)
    nullify(aux%aux_vars_bc)
  endif

#if 0
  if (associated(aux%aux_vars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call MiscibleAuxVarDestroy(aux%aux_vars_ss(iaux))
    enddo
    deallocate(aux%aux_vars_ss)
    nullify(aux%aux_vars_ss)
  endif
#endif

  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%miscible_parameter)) then
    if (associated(aux%miscible_parameter%dencpr)) deallocate(aux%miscible_parameter%dencpr)
    nullify(aux%miscible_parameter%dencpr)
    if (associated(aux%miscible_parameter%ckwet)) deallocate(aux%miscible_parameter%ckwet)
    nullify(aux%miscible_parameter%ckwet)
    if (associated(aux%miscible_parameter%sir)) deallocate(aux%miscible_parameter%sir)
    nullify(aux%miscible_parameter%sir)
    deallocate(aux%miscible_parameter)
  endif
  nullify(aux%miscible_parameter)
! if (associated(aux%res_old_AR)) deallocate(aux%res_old_AR)
! if (associated(aux%res_old_FL)) deallocate(aux%res_old_FL)
  if (associated(aux%delx)) deallocate(aux%delx)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine MiscibleAuxDestroy

end module Miscible_Aux_module
