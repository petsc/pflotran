 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! VERSION/REVISION HISTORY
 ! $Id: trdynmem.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
 ! $Log: trdynmem.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!
! Revision 1.3  2003/05/13 14:56:47  lichtner
! added header
!
!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 module trdynmem_module
  implicit none
  save
  public

#define PETSC_AVOID_DECLARATIONS
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscda.h"
#ifdef USE_PETSC216
    ! petsc-2.1.6
#include "include/finclude/petscsles.h"
#endif
#undef PETSC_AVOID_MPIF_H
#undef PETSC_AVOID_DECLARATIONS

!...........................................................
! connections, grid geometry
!...........................................................
    integer, pointer :: nd1(:), nd2(:)
    integer, pointer :: mblkbc(:), ibconn(:)
    real*8, pointer  :: dist1(:),dist2(:),distbc(:),area(:),areabc(:), &
                        vlbc(:), vgbc(:)! , xphibc(:)
    real*8, pointer  :: vb(:)

    integer, pointer :: nL2G(:), nG2L(:), nL2A(:)
    real*8, pointer :: ghost_loc_p(:)
    integer, allocatable :: nreg_val(:), nregbc(:)
    Vec :: nreg
!...........................................................
! Pitzer model
!...........................................................
  real*8, allocatable :: ah2o(:)

!...........................................................
! solver
!...........................................................
  Mat :: A
  Vec :: b, x
  real*8, pointer :: b_p(:)

!...........................................................
! grid
!...........................................................
  Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc
  Vec :: ghost, ghost_loc
!...........................................................
! primary species
!...........................................................
  Vec :: c, cc, ccloc, dc, ptran_c0
  real*8, allocatable :: psi(:,:), ppsi(:,:), dpsi(:,:,:)
  real*8, allocatable :: gam(:,:), cprev(:)
  real*8, allocatable :: rtoth2o(:), rtotco2(:),rrtoth2o(:), rrtotco2(:)
  real*8, pointer :: c_p(:), cloc_p(:)

!...........................................................
! secondary species
!...........................................................
  real*8, allocatable :: cx(:,:), gamx(:,:)

!...........................................................
! gas species
!...........................................................
  real*8, allocatable :: psig(:,:),ppsig(:,:),pgas(:,:), &
                         dpsig(:,:,:)

!...........................................................
! kinetic minerals
!...........................................................
  integer, allocatable  :: i1kin(:,:),i2kin(:,:), &
                           j1kin(:,:),j2kin(:,:), &
                           k1kin(:,:),k2kin(:,:)
  real*8, allocatable   :: phik_reg(:,:),surf_reg(:,:)

! real*8, allocatable   :: rkin(:,:), rrkin(:,:)
  real*8, pointer   :: rkin_p(:), rrkin_p(:)

  Vec :: phik, phik0, surf, surf0, rkin, rrkin
  real*8, pointer :: phik_p(:), phik0_p(:), surf_p(:), surf0_p(:)

!...........................................................
! surface complexation
!...........................................................
  real*8, allocatable :: siteden(:,:),csorpf(:,:),csorp(:,:), &
          ccsorp(:,:),alogpf(:,:)
  
!...........................................................
! ion exchange
!...........................................................
  real*8, allocatable :: cec(:,:),xex(:,:),xxex(:,:),glam(:,:)
  
!...........................................................
! porosity, saturation, tortuosity
!...........................................................
  Vec :: porosity, por, porloc, sat, sat_loc, ssat, ssat_loc, &
         tortuosity, tort_loc
  real*8, pointer :: por_p(:)

!...........................................................
! fields
!...........................................................
  real*8, allocatable :: rho(:), sionic(:)
  Vec     :: temp, temploc, press, pressloc

!...........................................................
! CO2 EOS
!...........................................................
  Vec :: xphi_co2, den_co2, xphi_co2_loc, den_co2_loc

!...........................................................
! velocity fields
!...........................................................
  real*8, allocatable :: vl(:), vg(:)

!...........................................................
! compression
!...........................................................
  
  integer, allocatable :: irow1(:),irow2(:),icol1(:),icol2(:)
  integer, allocatable :: dfill(:), ofill(:)        

  contains

    subroutine trallocate

    use ptran_global_module

    implicit none

    allocate(rho(ngmax))
    allocate(sionic(ngmax))
    if (iact == 6) allocate(ah2o(ngmax))
    allocate(psi(ncomp,ngmax))
    allocate(ppsi(ncomp,ngmax))
    allocate(dpsi(ncomp,ncomp,ngmax))
    allocate(gam(ncomp,ngmax))
    if (ngas > 0) allocate(pgas(ngas,ngmax))
!   allocate(sat(ngmax))
!   allocate(ssat(ngmax))
    if (iphase.eq.2 .or. iphase.eq.0) then
      allocate(psig(ncomp,ngmax))
      allocate(ppsig(ncomp,ngmax))
      allocate(dpsig(ncomp,ncomp,ngmax))
    endif
    if (ncmplx > 0) then
      allocate(cx(ncmplx,ngmax))
      allocate(gamx(ncmplx,ngmax))
    endif
    if (nexsolid > 0) then
      allocate(xex(nexmax,nlmax))
      allocate(xxex(nexmax,nlmax))
      allocate(glam(nexmax,nlmax))
    endif
    if (nsrfmin > 0) then
      allocate(siteden(nsrfmx,nlmax))
      allocate(csorpf(nsrfmx,nlmax))
      allocate(csorp(nsrfmx,nlmax))
      allocate(ccsorp(nsrfmx,nlmax))
    endif
!   if (nkin > 0) then
!     allocate(rkin(nkin,nlmax))
!     allocate(rrkin(nkin,nlmax))
!   endif
    allocate(rtoth2o(nlmax))
    allocate(rtotco2(nlmax))

    allocate(rrtoth2o(nlmax))
    allocate(rrtotco2(nlmax))
    
    !allocate(xphi_co2(nlmax))
    !allocate(den_co2(nlmax))
!   allocate(xphi_co2_bc())

!   allocate(alogpf(nscxmx,nlmax))

end subroutine trallocate
end module trdynmem_module
