module oil_eos_module
#include "finclude/petscsys.h"
  use petscsys

  private 

      
  PetscReal, parameter :: Rg=8.3145D3
  PetscReal, parameter :: den0= 863.983, p0= 17236893.25, vis0=3.0D-3, fmwoil= 142.D0
  
  public Vis_oil, oil_eos
  
 contains

! ************************************************************************** !

subroutine Vis_oil(p,t,viso, ierr)
  use PFLOTRAN_Constants_module

   implicit none
   
   PetscReal, intent(in) :: p,t
   PetscReal, intent(out):: viso
   PetscErrorCode :: ierr
   
  viso= vis0 
end subroutine Vis_oil

! ************************************************************************** !

subroutine PSAT_oil(t,psat,ierr)
  PetscReal t,psat
  PetscErrorCode :: ierr
end subroutine PSAT_oil

! ************************************************************************** !

subroutine oil_eos(t,p, x_mid, x_heavy, denoil, h, scale, ierr)
    implicit none

    PetscReal, intent(in) :: t   ! Temperature in centigrade
    PetscReal, intent(in) :: p   ! Pressure in Pascals
    PetscReal, intent(in) :: x_mid, x_heavy, scale
    PetscReal, intent(out) :: denoil
    PetscReal, intent(out) :: h
    PetscErrorCode, intent(out) :: ierr
    
!   PetscReal  Pc,Tc,w
!   PetscReal  a,b,alaph
!   PetscReal  Tk
    PetscReal  x_light

    ierr = 0
    x_light=1.D0- x_mid- x_heavy

! Tk=t+273.15
! Now only heavy oil formular runs, using data of n-Eicosane 
 !  x=1.D0  ! x=x_heavy
 !  Tc=775D0 !K 
 !  Pc=1.11D6 !Pa
 !  w= 0.71 !Acentric factor
 !  alaph= 1D0 + (1-(Tk/Tc)**0.5D0) * (0.37464 + w*(1.54226-0.26992*w))
 !  a = 0.45724* Rg*Rg*Tc*Tc/Pc *alaph
 !  b = 0.0778 * Rg *Tc /Pc

! p = 
  denoil =den0+ 2.24D-6*(p-p0)
  denoil = denoil/ fmwoil
  h =100.D0

end subroutine oil_eos


subroutine Meth_Eos
end subroutine Meth_Eos

  end module oil_eos_module
