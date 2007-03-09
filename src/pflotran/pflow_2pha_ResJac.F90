
! introduced grid variables: e_total :: 1 dof
!                            c_total :: grid%npricomp dof
!                            p_total :: 1 dof
!                            s_total :: (grid%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

#define PPRESSURE_LOC(n)   xx_loc_p(1+(n-1)*grid%ndof)
#define PPRESSURE(n)       xx_p(1+(n-1)*grid%ndof)
#define PRESSURE(n)        yy_p(1+(n-1)*grid%ndof)
#define TTEMP_LOC(n)       xx_loc_p(2+(n-1)*grid%ndof)
#define TTEMP(n)           xx_p(2+(n-1)*grid%ndof)
#define TEMP(n)            yy_p(2+(n-1)*grid%ndof) 
#define CCONC_LOC(n)       xx_loc_p(3+(n-1)*grid%ndof)
#define CCONC(n)           xx_p(3+(n-1)*grid%ndof)
#define CONC(n)            yy_p(3+(n-1)*grid%ndof)
#define SSATG_LOC(n)       xx_loc_p(4+(n-1)*grid%ndof)
#define SSATG(n)           xx_p(4+(n-1)*grid%ndof)
#define SATG(n)            yy_p(4+(n-1)*grid%ndof)


               
  module TTPHASE_module

  use pflow_gridtype_module

private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

  real*8, parameter :: formeps   = 5.D-6
  real*8, parameter :: eps       = 1.D-5
  real*8, parameter :: floweps   = 1.D-24
  real*8, parameter :: satcuteps = 1.D-5

  public TTPHASEResidual, TTPHASEJacobin, pflow_2phase_initaccum, &
  pflow_update_2phase,pflow_2phase_initadj,TTPHASESolutionBC, &
  TTPhase_Update_Reason,pflow_2phase_massbal

  contains

  subroutine TTPHASESolutionBC(snes,xx,grid)
    implicit none

    SNES, intent(in) :: snes
    Vec, intent(in) :: xx

    type(pflowGrid), intent(inout) :: grid
    integer ierr
    PetscScalar, pointer :: xx_p(:)
    integer n

    call VecGetArrayF90(xx, xx_p, ierr)
 
    do n = 1, grid%nlmax
       if(CCONC(n)>(1.D0 - 0.5D0*eps)) CCONC(n) = 1.D0 - eps
       if(CCONC(n)<(0.5D0*eps)) CCONC(n) =  eps
       if(SSATG(n)>(1.D0 - 0.5D0*eps)) SSATG(n) = 1.D0 - eps
       if(SSATG(n)<(0.5D0*eps)) SSATG(n) =  eps
    end do

    call VecRestoreArrayF90(xx, xx_p, ierr)

  end subroutine TTPHASESolutionBC


  subroutine TTPhase_Update(xx,grid)
  
  use water_eos_module
  
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  Vec, intent(in) :: xx

  PetscScalar, pointer :: xx_p(:), yy_p(:), iphase_p(:)
  integer :: n,n0
  integer :: ierr,iipha

  real*8 :: pco2,p,tmp,xw
  
#include "definitions.h"

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)

  do n = 1,grid%nlmax
    n0=(n-1)* grid%ndof
    iipha=iphase_p(n)
        
    p = xx_p(n0+1)
       !call Henry_co2_noderiv( xx_p(n0+2),ppsat,ierr)
   !    pco2_sat= p * xx
    pco2 = p * xx_p(n0+3)
    xw = 1.D0 - xx_p(n0+3)

    if (xw > grid%yh2o_in_co2 .and. iipha==4 )then
      write(*,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      iphase_p(n) = 6
      xx_p(n0+4)=1.D0- formeps
      xx_p(n0+3)=1.D0 
        !  xx_p(n0+2)= yy_p(n0+2)
    end if
       ! gas ->  2ph 

    if (xx_p(n0+3) > 1.025D0  .and. iipha==2) then
      write(*,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      write(IUNIT2,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      iphase_p(n) = 6
      tmp= formeps
           ! if(tmp>1.D-6)tmp=1.D-6
     !     if(tmp<eps) tmp=eps
      xx_p(n0+4)=tmp
      tmp=1.D0 - grid%yh2o_in_co2 !1.D0+( xx_p(n0+3) -1.D0)*0.5D0 &
      !xx_p(n0+3)
          ! if(tmp>1.01D0)  tmp=1.01D0 !-0.75D0*ppsat/p
          ! if(xx_p(n0+3)<tmp) 
      xx_p(n0+3)=tmp
    end if
      
    if(xx_p(n0 + 4) >= (1.D0 - eps).and. iipha==6 )then
       ! print *,'S upcorin ',n,xx_p(n0+1:n0+4)
       
       ! xx_p(n0 + 3) = 1.D0 - 
      xx_p(n0 + 4) = 1.00D0  
      if(xw<eps)then
        iphase_p(n) = 4 ! 2ph ->gas 
        tmp = 1.D0 !- ppsat / p
        xx_p(n0 + 4) = 1.00D0
        !if(xx_p(n0 + 3)>1.D0)
        xx_p(n0 + 3) =1.D0
       ! if(xx_p(n0 + 3)< 0.D0)xx_p(n0 + 3)=0.D0
        write(*,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
        write(IUNIT2,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      end if
    endif

    if(xx_p(n0 + 4)<= formeps .and. iipha==6 )then
       ! print *,'S lowcorin ',n,xx_p(n0+1:n0+4)
   
      !if( xx_p(n0+3)< 1.D0)then
      write(*,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      write(IUNIT2,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+4)
      iphase_p(n) = 2 ! 2ph -> Liq
      tmp = 1.D0 !- grid%yh2o_in_co2 !.95D0 - ppsat / p
        !if(  xx_p(n0 + 3) > tmp)
      xx_p(n0 + 3)=tmp
      xx_p(n0 + 4) = 0.D0
      !end if
    !  if( xx_p(n0 + 4) <0.D0 ) xx_p(n0 + 4) =  0.D0
    endif
      
   !  if(xx_p(n0 + 3) < 0.D0) xx_p(n0 + 3) = 0.D0
   !   if(xx_p(n0 + 3) > 1.D0) xx_p(n0 + 3) =  yy_p(n0 + 3)

  end do

  !print *,iphase_p
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)

  end subroutine TTPhase_Update


  subroutine TTPhase_Update_Reason(reason,grid)
  
  implicit none
 
  integer, intent(out):: reason
  type(pflowGrid), intent(inout) :: grid
  PetscScalar, pointer :: xx_p(:)
  integer :: n,n0,re
  integer re0, ierr

  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)

  re=1
  do n = 1,grid%nlmax
     n0=(n-1)* grid%ndof
     if(xx_p(n0 + 3) > 1.05D0) re=0; exit
     if(xx_p(n0 + 3) < -.05D0) re=0; exit
     if(xx_p(n0 + 4) > 1.01D0) re=0; exit
     if(xx_p(n0 + 4) < -1.D-2) re=0; exit
  end do
  call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  
  call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  
  if(grid%commsize >1)then
    call MPI_REDUCE(re, re0,1, MPI_INTEGER,MPI_SUM,0, &
    PETSC_COMM_WORLD,ierr)
    call MPI_BCAST(re0,1, MPI_INTEGER, 0,PETSC_COMM_WORLD,ierr)
    if(re0<grid%commsize) re=0
  endif
  reason=re
  
  if(reason<=0) print *,'Sat or Con out of Region'
  
  end subroutine TTPhase_Update_Reason




  subroutine TTPHASEResidual(snes,xx,r,grid,ierr)

    use water_eos_module
    use co2eos_module
    use mixture_module
    use span_wagner_module

    implicit none
 
    SNES, intent(in) :: snes
    Vec, intent(inout) :: xx
    Vec, intent(out) :: r
    type(pflowGrid), intent(inout) :: grid

 
  integer :: ierr
  integer :: n, ng, nc, nr
  integer :: i, i1, i2, j, jn, jng, jm1, jm2, jmu
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2, p1, p2, t1, t2, c1, c2,&
             s1, s2
  integer :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
  integer :: i1_hencoeff, i2_hencoeff
  integer :: ibc  ! Index that specifies a boundary condition block
  
! real*8 :: term1, term2, term3


  PetscScalar, pointer ::accum_p(:)

  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), &
               ddensity_p(:), ddensity_loc_p(:),&
               phis_p(:),  &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), &
               d_p_p(:), d_p_loc_p(:), &
               d_t_p(:), d_t_loc_p(:), &
               d_c_p(:), d_c_loc_p(:), &
               d_s_p(:), d_s_loc_p(:), &
               avgmw_p(:),avgmw_loc_p(:),avgmw_c_p(:),avgmw_c_loc_p(:),&
               hh_p(:), hh_loc_p(:), &
               h_p_p(:), h_p_loc_p(:), &
               h_t_p(:), h_t_loc_p(:), &
               h_c_p(:), h_c_loc_p(:), &
               h_s_p(:), h_s_loc_p(:), &
               uu_p(:),       &
               u_p_p(:), u_t_p(:), u_c_p(:), u_s_p(:), &
               hen_p(:),  hen_loc_p(:), &  
               hen_p_p(:), hen_p_loc_p(:), &
               hen_t_p(:), hen_t_loc_p(:), &
               hen_c_p(:), hen_c_loc_p(:), &
               hen_s_p(:), hen_s_loc_p(:), &
               df_p(:), df_loc_p(:), &
               df_p_p(:), df_p_loc_p(:),&
               df_t_p(:), df_t_loc_p(:),&
               df_c_p(:), df_c_loc_p(:),&
               df_s_p(:), df_s_loc_p(:)
               
  PetscScalar, pointer :: pc_p(:), pc_loc_p(:),&
                          pc_p_p(:), pc_p_loc_p(:),&
                          pc_t_p(:), pc_t_loc_p(:),&
                          pc_c_p(:), pc_c_loc_p(:),&
                          pc_s_p(:), pc_s_loc_p(:),&
                          kvr_p(:), kvr_loc_p(:),&
                          kvr_p_p(:), kvr_p_loc_p(:),&
                          kvr_t_p(:), kvr_t_loc_p(:),&
                          kvr_c_p(:), kvr_c_loc_p(:),&
                          kvr_s_p(:), kvr_s_loc_p(:)

  PetscScalar, pointer :: iphase_loc_p(:),icap_p(:),&
                          icap_loc_p(:), ithrm_loc_p(:)

  integer :: iicap,iiphase

  real*8 :: dd1, dd2, eng, cond, den, &
            eengl,eengg, &
            fluxcl,fluxcg,fluxe, fluxh, flux, gravity, fluxl,&
            fluxlh,fluxlv, fluxg,fluxgh,fluxgv, fluxv, q,  &
            v_darcy,hflx,pvoldt, voldt, accum, pvol
  real*8 :: p_vapor,qu_rate,SSATW
  real*8 :: dd, f1, f2, ff, por1, por2, perm1, perm2
  real*8 :: Dphi,D0
  real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants at upstream, downstream faces.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol,density_ave,difff,diffl,diffg,difluxl,difluxg
  real*8 :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 !, qqsrc
  real*8 :: cw,cw1,cw2, xxlw,xxla,xxgw,xxga
  real*8 :: upweight
  real*8 :: ukvr,uhh,uconc
  real*8 :: dddt,dddp,fg,dfgdp,dfgdt,dhdt,dhdp,dvdt,dvdp, rho, visc
 
  grid%vvlbc=0.D0
  grid%vvgbc=0.D0
  grid%vvl_loc=0.D0
  grid%vvg_loc=0.D0

 !call  TTPhase_Update(xx,grid)
  
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecGetArrayF90(grid%hh, hh_p, ierr)
  call VecGetArrayF90(grid%uu, uu_p, ierr)
  call VecGetArrayF90(grid%df, df_p, ierr)
  call VecGetArrayF90(grid%hen, hen_p,ierr) 
  call VecGetArrayF90(grid%pcw, pc_p, ierr)
  call VecGetArrayF90(grid%kvr, kvr_p, ierr)
  call VecGetArrayF90(grid%icap,icap_p,ierr)
! call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)
! --- calculate derivative when evaluating residual---
  if (grid%ideriv == 1) then
    call VecGetArrayF90(grid%d_p, d_p_p, ierr)
    call VecGetArrayF90(grid%d_t, d_t_p, ierr)
    call VecGetArrayF90(grid%d_c, d_c_p, ierr)
    call VecGetArrayF90(grid%d_s, d_s_p, ierr)
    call VecGetArrayF90(grid%avgmw_c,avgmw_c_p, ierr)
    call VecGetArrayF90(grid%h_p, h_p_p, ierr)
    call VecGetArrayF90(grid%h_t, h_t_p, ierr)
    call VecGetArrayF90(grid%h_c, h_c_p, ierr)
    call VecGetArrayF90(grid%h_s, h_s_p, ierr)
    call VecGetArrayF90(grid%u_p, u_p_p, ierr)
    call VecGetArrayF90(grid%u_t, u_t_p, ierr)
    call VecGetArrayF90(grid%u_c, u_c_p, ierr)
    call VecGetArrayF90(grid%u_s, u_s_p, ierr)
    call VecGetArrayF90(grid%df_p, df_p_p, ierr)
    call VecGetArrayF90(grid%df_t, df_t_p, ierr)
    call VecGetArrayF90(grid%df_c, df_c_p, ierr)
    call VecGetArrayF90(grid%df_s, df_s_p, ierr)
    call VecGetArrayF90(grid%hen_p, hen_p_p, ierr)
    call VecGetArrayF90(grid%hen_t, hen_t_p, ierr)
    call VecGetArrayF90(grid%hen_c, hen_c_p, ierr)
    call VecGetArrayF90(grid%hen_s, hen_s_p, ierr)
    call VecGetArrayF90(grid%pc_p, pc_p_p, ierr)
    call VecGetArrayF90(grid%pc_t, pc_t_p, ierr)
    call VecGetArrayF90(grid%pc_c, pc_c_p, ierr)
    call VecGetArrayF90(grid%pc_s, pc_s_p, ierr)
    call VecGetArrayF90(grid%kvr_p, kvr_p_p, ierr)
    call VecGetArrayF90(grid%kvr_t, kvr_t_p, ierr)
    call VecGetArrayF90(grid%kvr_c, kvr_c_p, ierr)
    call VecGetArrayF90(grid%kvr_s, kvr_s_p, ierr)
  endif
!------------------------------------------------------ 

!-----  phase properities ---- last time step---
  do n = 1, grid%nlmax
    jn = 1 + (n-1)*grid%nphase
       
    ii1=jn  !1+(n-1)*grid%nphase
    ii2=n*grid%nphase
    iicap=icap_p(n)

    if (grid%ideriv .eq. 1) then

      call mixture_eos(PPRESSURE(n),TTEMP(n),CCONC(n),SSATG(n),&
        grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
        iicap, grid%swir(iicap),&
        grid%lambda(iicap),grid%alpha(iicap),&
        grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
        ddensity_p(ii1:ii2),d_p_p(ii1:ii2),d_t_p(ii1:ii2),&
        d_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        d_s_p(ii1:ii2),avgmw_p(ii1:ii2),&
        avgmw_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        hh_p(ii1:ii2),h_p_p(ii1:ii2),h_t_p(ii1:ii2),&
        h_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        h_s_p(ii1:ii2),&
        uu_p(ii1:ii2),u_p_p(ii1:ii2),u_t_p(ii1:ii2),&
        u_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        u_s_p(ii1:ii2),&
        df_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase *grid%nspec),&
        df_p_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        df_t_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        df_c_p(1+(n-1)*grid%nphase*grid%nspec*grid%npricomp:&
               n*grid%nphase*grid%nspec*grid%npricomp),&
        df_s_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_p_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_t_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_c_p(1+(n-1)*grid%nphase*grid%nspec*grid%npricomp:&
                n*grid%nphase*grid%nspec*grid%npricomp),&
        hen_s_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        grid%xxphi_co2(n),&
        pc_p(ii1:ii2),pc_p_p(ii1:ii2),pc_t_p(ii1:ii2),&
        pc_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        pc_s_p(ii1:ii2),& 
        kvr_p(ii1:ii2),kvr_p_p(ii1:ii2),kvr_t_p(ii1:ii2),&
        kvr_c_p(1+(n-1)*grid%nphase*grid%npricomp:n*grid%nphase*grid%npricomp),&
        kvr_s_p(ii1:ii2),&
        ierr,grid%itable)
    else     
      call mixture_eos_noderiv (PPRESSURE(n),TTEMP(n),CCONC(n),SSATG(n),&
        grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
        iicap, grid%swir(iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
        grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
        ddensity_p(ii1:ii2),avgmw_p(ii1:112),hh_p(ii1:ii2),uu_p(ii1:ii2),&
        df_p((1+(n-1)*grid%nphase*grid%nspec):(n*grid%nphase*grid%nspec)),&
        hen_p((1+(n-1)*grid%nphase*grid%nspec):(n*grid%nphase*grid%nspec)),&    
        grid%xxphi_co2(n),pc_p(ii1:ii2),kvr_p(ii1:ii2),&
        ierr,grid%itable)
    endif
    
!   if(n < 5) print *,'pflow_2ph: ',n,grid%ideriv,grid%xxphi_co2(n),cconc(n)
  enddo

  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%ddensity, ddensity_p, ierr)

  call VecRestoreArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecRestoreArrayF90(grid%hh, hh_p, ierr)
  call VecRestoreArrayF90(grid%uu, uu_p, ierr)
  call VecRestoreArrayF90(grid%df, df_p, ierr)
  call VecRestoreArrayF90(grid%hen, hen_p,ierr)
  call VecRestoreArrayF90(grid%pcw,pc_p,ierr)
  call VecRestoreArrayF90(grid%kvr, kvr_p, ierr)
! call VecRestoreArrayF90(grid%iphase,iphase_p,ierr)
  
  call VecRestoreArrayF90(grid%icap,icap_p,ierr)

  if (grid%ideriv == 1) then  
    call VecRestoreArrayF90(grid%d_p, d_p_p, ierr)
    call VecRestoreArrayF90(grid%d_t, d_t_p, ierr)
    call VecRestoreArrayF90(grid%d_c, d_c_p, ierr)
    call VecRestoreArrayF90(grid%d_s, d_s_p, ierr)
    call VecRestoreArrayF90(grid%avgmw_c,avgmw_c_p, ierr)
    call VecRestoreArrayF90(grid%h_p, h_p_p, ierr)
    call VecRestoreArrayF90(grid%h_t, h_t_p, ierr)
    call VecRestoreArrayF90(grid%h_c, h_c_p, ierr)
    call VecRestoreArrayF90(grid%h_s, h_s_p, ierr)
    call VecRestoreArrayF90(grid%u_p, u_p_p, ierr)
    call VecRestoreArrayF90(grid%u_t, u_t_p, ierr)
    call VecRestoreArrayF90(grid%u_c, u_c_p, ierr)
    call VecRestoreArrayF90(grid%u_s, u_s_p, ierr)
    call VecRestoreArrayF90(grid%df_p, df_p_p, ierr)
    call VecRestoreArrayF90(grid%df_t, df_t_p, ierr)
    call VecRestoreArrayF90(grid%df_c, df_c_p, ierr)
    call VecRestoreArrayF90(grid%df_s, df_s_p, ierr)
    call VecRestoreArrayF90(grid%hen_p, hen_p_p, ierr)
    call VecRestoreArrayF90(grid%hen_t, hen_t_p, ierr)
    call VecRestoreArrayF90(grid%hen_c, hen_c_p, ierr)
    call VecRestoreArrayF90(grid%hen_s, hen_s_p, ierr)
    call VecRestoreArrayF90(grid%pc_p, pc_p_p, ierr)
    call VecRestoreArrayF90(grid%pc_t, pc_t_p, ierr)
    call VecRestoreArrayF90(grid%pc_c, pc_c_p, ierr)
    call VecRestoreArrayF90(grid%pc_s, pc_s_p, ierr)
    call VecRestoreArrayF90(grid%kvr_p, kvr_p_p, ierr)
    call VecRestoreArrayF90(grid%kvr_t, kvr_t_p, ierr)
    call VecRestoreArrayF90(grid%kvr_c, kvr_c_p, ierr)
    call VecRestoreArrayF90(grid%kvr_s, kvr_s_p, ierr)
  endif


  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%ddensity, &
                            INSERT_VALUES,grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof,grid%ddensity,INSERT_VALUES, &
                          grid%ddensity_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%avgmw,INSERT_VALUES, &
                            grid%avgmw_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%avgmw, INSERT_VALUES, &
                          grid%avgmw_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                            grid%hh_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                          grid%hh_loc, ierr)
!  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%uu, INSERT_VALUES, &
 !                           grid%uu_loc, ierr)
 ! call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%uu, INSERT_VALUES, &
 !                         grid%uu_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_NphaNspec_dof,grid%df,INSERT_VALUES, &
                            grid%df_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_NphaNspec_dof, grid%df, INSERT_VALUES, &
                          grid%df_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_NphaNspec_dof,grid%hen,INSERT_VALUES, &
                            grid%hen_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_NphaNspec_dof, grid%hen, INSERT_VALUES, &
                          grid%hen_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%pcw, INSERT_VALUES, &
                            grid%pcw_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%pcw, INSERT_VALUES, &
                          grid%pcw_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%kvr, INSERT_VALUES, &
                            grid%kvr_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%kvr, INSERT_VALUES, &
                          grid%kvr_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof,grid%porosity,INSERT_VALUES, &
                            grid%porosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                          grid%porosity_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_xx, &
                            INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_xx, &
                          INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_yy, &
                            INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_yy, &
                          INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_zz, &
                            INSERT_VALUES, grid%perm_zz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_zz, &
                          INSERT_VALUES, grid%perm_zz_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%ithrm, &
                            INSERT_VALUES, grid%ithrm_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%ithrm, &
                          INSERT_VALUES, grid%ithrm_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%iphas, &
                           INSERT_VALUES, grid%iphas_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%iphas, &
                          INSERT_VALUES, grid%iphas_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%icap, &
                            INSERT_VALUES, grid%icap_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%icap, &
                          INSERT_VALUES, grid%icap_loc, ierr)

! now assign access pointer to local variables
  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
! call VecGetArrayF90(grid%yy, yy_p, ierr)
 
 
  ! notice:: here we assume porosity is constant

  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%avgmw_loc, avgmw_loc_p, ierr)
! call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(grid%uu, uu_p, ierr)
! call VecGetArrayF90(grid%h, h_p, ierr)

  call VecGetArrayF90(grid%df_loc, df_loc_p, ierr)
  call VecGetArrayF90(grid%hen_loc, hen_loc_p, ierr)
  call VecGetArrayF90(grid%pcw_loc, pc_loc_p, ierr)
  call VecGetArrayF90(grid%kvr_loc, kvr_loc_p, ierr)

  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(grid%vl, vl_p, ierr)
  call VecGetArrayF90(grid%iphas_loc, iphase_loc_p, ierr)
 !print *,' Finished scattering non deriv'


  if (grid%ideriv == 1) then
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%d_p,INSERT_VALUES, &
                            grid%d_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                          grid%d_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%d_t,INSERT_VALUES, &
                            grid%d_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_t, INSERT_VALUES, &
                          grid%d_t_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphancomp_dof,grid%d_c, &
                            INSERT_VALUES, grid%d_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphancomp_dof,grid%d_c,INSERT_VALUES, &
                          grid%d_c_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%d_s,INSERT_VALUES, &
                            grid%d_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_s, INSERT_VALUES, &
                          grid%d_s_loc, ierr)
!print *,' Finished scattering deriv: d'

    call DAGlobalToLocalBegin(grid%da_nphancomp_dof, grid%avgmw_c, &
                            INSERT_VALUES, grid%avgmw_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphancomp_dof, grid%avgmw_c, &
                          INSERT_VALUES, grid%avgmw_c_loc, ierr)
  !print *,' Finished scattering deriv: a'
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%h_p,INSERT_VALUES, &
                            grid%h_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_p, INSERT_VALUES, &
                          grid%h_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%h_t,INSERT_VALUES, &
                            grid%h_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_t, INSERT_VALUES, &
                          grid%h_t_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphancomp_dof,grid%h_c, &
                            INSERT_VALUES,grid%h_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphancomp_dof,grid%h_c,INSERT_VALUES, &
                          grid%h_c_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%h_s,INSERT_VALUES, &
                            grid%h_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_s, INSERT_VALUES, &
                          grid%h_s_loc, ierr)

    call DAGlobalToLocalBegin(grid%da_nphanspec_dof,grid%df_p, &
                            INSERT_VALUES,grid%df_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof,grid%df_p, &
                          INSERT_VALUES,grid%df_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspec_dof,grid%df_t, &
                            INSERT_VALUES,grid%df_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof, grid%df_t, &
                          INSERT_VALUES, grid%df_t_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspec_dof, grid%df_s, &
                            INSERT_VALUES, grid%df_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof, grid%df_s, &
                          INSERT_VALUES, grid%df_s_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspecncomp_dof, grid%df_c, &
                            INSERT_VALUES, grid%df_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspecncomp_dof, grid%df_c, &
                          INSERT_VALUES, grid%df_c_loc, ierr)

    call DAGlobalToLocalBegin(grid%da_nphanspec_dof, grid%hen_p, &
                            INSERT_VALUES, grid%hen_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof, grid%hen_p, &
                          INSERT_VALUES, grid%hen_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspec_dof, grid%hen_t, &
                            INSERT_VALUES, grid%hen_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof, grid%hen_t, &
                          INSERT_VALUES, grid%hen_t_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspec_dof, grid%hen_s, &
                            INSERT_VALUES, grid%hen_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspec_dof, grid%hen_s, &
                          INSERT_VALUES, grid%hen_s_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphanspecncomp_dof, grid%hen_c, &
                            INSERT_VALUES, grid%hen_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphanspecncomp_dof, grid%hen_c, &
                          INSERT_VALUES, grid%hen_c_loc, ierr)
!print *,' Finished scattering deriv: df, hen'
   
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%pc_p,INSERT_VALUES, &
                            grid%pc_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%pc_p, INSERT_VALUES, &
                          grid%pc_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%pc_t,INSERT_VALUES, &
                            grid%pc_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%pc_t, INSERT_VALUES, &
                          grid%pc_t_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphancomp_dof, grid%pc_c, &
                            INSERT_VALUES, grid%pc_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphancomp_dof, grid%pc_c, &
                          INSERT_VALUES, grid%pc_c_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%pc_s,INSERT_VALUES, &
                            grid%pc_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%pc_s, INSERT_VALUES, &
                          grid%pc_s_loc, ierr)


    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%kvr_p, &
                            INSERT_VALUES, grid%kvr_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof,grid%kvr_p,INSERT_VALUES, &
                          grid%kvr_p_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%kvr_t, &
                            INSERT_VALUES, grid%kvr_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof,grid%kvr_t,INSERT_VALUES, &
                          grid%kvr_t_loc, ierr)  

    call DAGlobalToLocalBegin(grid%da_nphancomp_dof, grid%kvr_c, &
                            INSERT_VALUES, grid%kvr_c_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphancomp_dof, grid%kvr_c, &
                          INSERT_VALUES, grid%kvr_c_loc, ierr)
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%kvr_s, &
                            INSERT_VALUES, grid%kvr_s_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof,grid%kvr_s,INSERT_VALUES, &
                          grid%kvr_s_loc, ierr)

    call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecGetArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecGetArrayF90(grid%d_c_loc, d_c_loc_p, ierr)
    call VecGetArrayF90(grid%d_s_loc, d_s_loc_p, ierr)
    call VecGetArrayF90(grid%avgmw_c_loc,avgmw_c_loc_p, ierr)
    call VecGetArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecGetArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecGetArrayF90(grid%h_s_loc, h_c_loc_p, ierr)
    call VecGetArrayF90(grid%h_s_loc, h_s_loc_p, ierr)
    call VecGetArrayF90(grid%df_p_loc, df_p_loc_p, ierr)
    call VecGetArrayF90(grid%df_t_loc, df_t_loc_p, ierr)
    call VecGetArrayF90(grid%df_c_loc, df_c_loc_p, ierr)
    call VecGetArrayF90(grid%df_s_loc, df_s_loc_p, ierr)
    call VecGetArrayF90(grid%hen_p_loc, hen_p_loc_p, ierr)
    call VecGetArrayF90(grid%hen_t_loc, hen_t_loc_p, ierr)
    call VecGetArrayF90(grid%hen_c_loc, hen_c_loc_p, ierr)
    call VecGetArrayF90(grid%hen_s_loc, hen_s_loc_p, ierr)
    call VecGetArrayF90(grid%pc_p_loc, pc_p_loc_p, ierr)
    call VecGetArrayF90(grid%pc_t_loc, pc_t_loc_p, ierr)
    call VecGetArrayF90(grid%pc_c_loc, pc_c_loc_p, ierr)
    call VecGetArrayF90(grid%pc_s_loc, pc_s_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_p_loc, kvr_p_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_t_loc, kvr_t_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_c_loc, kvr_c_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_s_loc, kvr_s_loc_p, ierr)

  endif

  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
  endif

 !accum_p =0.D0
 !print *,' Finished scattering  deriv'
  r_p=0.D0
!--------------------------------------------------------------------------
! Calculate accumulation term for interior and exterior nodes.
!--------------------------------------------------------------------------
! print *,grid%rtot
  
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1
    s1 = c1 + grid%npricomp

    pvol = volume_p(n)*porosity_loc_p(ng)
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt
    SSATW = 1.D0-SSATG_LOC(ng)
    iiphase = iphase_loc_p(ng)
    
!   Pressure equation accumulation term
    accum = 0.d0
!   do j = 1, grid%nphase
      j=1
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
      
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
      r_p(p1)=-accum_p(p1)/grid%dt  

!     units: volume_p [m^3]; ddensity [mol/dm^3]; r [kmol/s]
      if(iiphase==6) then ! have 2ph
        r_p(p1)=r_p(p1) + pvoldt * ddensity_loc_p(1+jng) * SSATG_LOC(ng) ! kmol/s
        r_p(p1)=r_p(p1) + pvoldt * ddensity_loc_p(jng)*SSATW
        r_p(p1)=r_p(p1) - grid%rtot(n,1) - grid%rtot(n,2)
!       print *,'build r liq',n,ng,volume_p(n),ddensity_loc_p(jng), &
!       ddensity_loc_p(1+jng),SSATG_LOC(ng),SSATW,porosity_loc_p(ng),r_p(p1)
      end if

      if(iiphase==2) then ! have liq
        r_p(p1)=r_p(p1) + pvoldt *ddensity_loc_p(jng) &
        - grid%rtot(n,1) - grid%rtot(n,2)
      endif
      if(iiphase==4) then ! have gas
        r_p(p1)=r_p(p1) + pvoldt * ddensity_loc_p(1+jng)
      endif
!   enddo

!   Heat equation accumulation term
    i = ithrm_loc_p(ng)
!   j = grid%jh2o
    j=1
    jn = j+(n-1)*grid%nphase
    jng = j+(ng-1)*grid%nphase
    ! rho U = rho H - p, and here the eng has unit as J*m^(-3)
    eengl = ddensity_loc_p(jng)*uu_p(jn)
    eengg = ddensity_loc_p(1+jng)*uu_p(1+jn)
 !  eengl=eengl * SSATW
 !  eengg=eengg * SSATG_LOC(ng)
 
 
    r_p(t1) = ((1.d0-porosity_loc_p(ng)) * grid%dencpr(i) * &
                TTEMP_LOC(ng) * volume_p(n) - accum_p(t1))/grid%dt 
    if(iiphase==6) then ! have 2ph
      r_p(t1) = r_p(t1) + pvoldt *( eengg * SSATG_LOC(ng)+ eengl * SSATW)
    end if

    if(iiphase==4) then ! have gas
      r_p(t1) = r_p(t1) + pvoldt * eengg 
    end if
    if(iiphase==2) then ! have liq
      r_p(t1)=r_p(t1) +  pvoldt * eengl 
    end if

   !call Henry_coeff(PPRESSURE_LOC(ng),TTEMP_LOC(ng),henrycoeff)
    xxga = CCONC_LOC(ng)
    xxgw = 1.D0 - xxga
    xxla = xxga*Hen_loc_p(2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec)
    xxlw = 1.D0 - xxla

    r_p(s1)= -accum_p(s1)/grid%dt
    if(iiphase==6) then ! have gas
      r_p(s1) = r_p(s1) + pvoldt * (ddensity_loc_p(1 + jng) * xxgw  *  &
           SSATG_LOC(ng) + ddensity_loc_p(jng) * xxlw * SSATW)
      r_p(s1) = r_p(s1) - grid%rtot(n,1)
    end if

    if(iiphase==4) then ! have gas
      r_p(s1) = r_p(s1)+ pvoldt * ddensity_loc_p(1 + jng) * xxgw 
    end if
    if(iiphase==2) then ! have liq
      r_p(s1) = r_p(s1)+ pvoldt* ddensity_loc_p(jng) * xxlw
      r_p(s1) = r_p(s1) - grid%rtot(n,1)
    ! Liq phase water component mass balance
      !print *,'Res S ',n, xxlw ,xxla,  ddensity_loc_p(jng),pvoldt,&
      !      Hen_loc_p(2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec)
    endif   
 


    if(iiphase==6)then
      accum=accum_p(c1)
      r_p(c1) = (pvol * ddensity_loc_p(1+jng) * xxgw * SSATG_LOC(ng) &
           * grid%ret-accum)/grid%dt
     ! Gas phase water component mass balance  
    else
      if(iiphase==2) r_p(c1) = SSATG_LOC(ng)
      if(iiphase==4) r_p(c1) = SSATG_LOC(ng)-1.D0
    endif


!   Reaction term contribution to concentration equation
    qu_rate=0.D0  ! evaporation get positive
    if(iphase_loc_p(ng)==6) then
      p_vapor=PPRESSURE_LOC(ng)*xxgw
    
      call psat(TTEMP_LOC(ng),sat_pressure,ierr)

      if( SSATG_LOC(ng)<(1.D0-eps).and. SSATG_LOC(ng)>eps )then
        qu_rate=-grid%qu_kin*(p_vapor- &
        grid%yh2o_in_co2*PPRESSURE_LOC(ng))!*xxlw)!*&
      end if
      r_p(c1)=r_p(c1) - qu_rate
    end if
  end do

!************************************************************************
 ! add source/sink terms
 
  do nr = 1, grid%nblksrc
      
    kk1 = grid%k1src(nr) - grid%nzs
    kk2 = grid%k2src(nr) - grid%nzs
    jj1 = grid%j1src(nr) - grid%nys
    jj2 = grid%j2src(nr) - grid%nys
    ii1 = grid%i1src(nr) - grid%nxs
    ii2 = grid%i2src(nr) - grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1 = grid%qsrc(i,nr)
        csrc1 = grid%csrc(i,nr)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1 = f1*grid%qsrc(i,nr) + f2*grid%qsrc(i-1,nr)
        csrc1 = f1*grid%csrc(i,nr) + f2*grid%csrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
   !print *,'pflow2ph : ', grid%myrank,i,grid%timesrc(i,nr), &
   !grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1,tsrc1
 
    qsrc1 = qsrc1 / grid%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / grid%fmwco2
  
  ! Here assuming regular mixture injection. i.e. no extra H from mixing 
  ! within injected fluid.
  
    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)
            p1 = 1+(n-1)*grid%ndof
            t1 = p1 + 1
            c1 = t1 + 1
            s1 = c1 + 1

            call wateos_noderiv(tsrc1,PPRESSURE_LOC(ng),dw_kg,dw_mol, &
            enth_src_h2o,grid%scale,ierr)

!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]

!           qqsrc = qsrc1/dw_mol ! [kmol/s / mol/dm^3 = kmol/m^3]
              
            r_p(p1) = r_p(p1) - qsrc1 
            r_p(t1) = r_p(t1) - qsrc1*enth_src_h2o
            r_p(s1) = r_p(s1) - qsrc1

!           print *,'pflow2ph_h2o: ',nr,n,ng,tsrc1,dw_mol,dw_mol*grid%fmwh2o, &
!           qsrc1
          enddo
        enddo
      enddo
    endif  
    
    if (csrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
            n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
            ng = grid%nL2G(n)
            jng= 2 + (ng-1)*grid%nphase
            p1 = 1+(n-1)*grid%ndof
            t1 = p1 + 1
            !c1 = t1 + 1
        
!           duan eos
!           call duanco2(tsrc1,PPRESSURE_LOC(ng)/1D5,dco2,fugco2,co2_phi)
!           call ENTHALPY(tsrc1+273.15D0,1.D-3/dco2,1.D0/co2_phi, &
!           enth_src_co2)
!           enth_src_co2=enth_src_co2 * 1.D-3     
 
         !  span-wagner
            rho = ddensity_loc_p(jng)*grid%fmwco2 
            call co2_span_wagner(PPRESSURE_LOC(ng)*1.D-6,tsrc1+273.15D0, &
            rho,dddt,dddp,fg,dfgdp,dfgdt, &
            eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,grid%itable)

         !  units: rho [kg/m^3]; csrc1 [kmol/s]

            enth_src_co2 = enth_src_co2 * grid%fmwco2

            r_p(p1) = r_p(p1) - csrc1
            r_p(t1) = r_p(t1) - csrc1 * enth_src_co2
           !r_p(s1) = r_p(s1) - csrc1

           print *,'pflow2ph_co2: ',nr,n,ng,tsrc1,rho,grid%fmwco2,csrc1
          enddo
        enddo
      enddo
    endif
  
  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
      
    !  do kk = kk1, kk2
     !   do jj = jj1, jj2
       !   do ii = ii1, ii2
          !    n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
           !   ng = grid%nL2G(n)
           !   p1 = 1+(n-1)*grid%ndof
            !  t1 = p1 + 1
           !   c1 = t1 + 1
          !    qqsrc = qsrc1/ddensity_loc_p(ng)
         !     enth_src = hh_loc_p(ng)
        !      r_p(p1) = r_p(p1) - qsrc1
       !       r_p(t1) = r_p(t1) - qsrc1*enth_src
      !        r_p(c1) = r_p(c1) - qqsrc*CCONC_LOC(ng)
     !     enddo
    !    enddo
   !   enddo
  !  endif
  enddo


!*********************************************************************


 
! stop
!---------------------------------------------------------------------------
! Flux terms for interior nodes
! Be careful here, we have velocity field for every phase
!---------------------------------------------------------------------------
 
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   

    p1 = 1 + (n1-1)*grid%ndof; t1 = p1+1; c1 = t1+1; s1=c1 + grid%npricomp
    p2 = 1 + (n2-1)*grid%ndof; t2 = p2+1; c2 = t2+1; s2=c2 + grid%npricomp
   
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)


    select case(ip1)
    case(1) 
       perm1 = perm_xx_loc_p(m1)
    case(2)
       perm1 = perm_yy_loc_p(m1)
    case(3)
       perm1 = perm_zz_loc_p(m1)
    end select
    
    select case(ip2)
    case(1) 
       perm2 = perm_xx_loc_p(m2)
    case(2)
       perm2 = perm_yy_loc_p(m2)
    case(3)
       perm2 = perm_zz_loc_p(m2)
    end select


    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    fluxl = 0.d0
    fluxlh = 0.d0
    fluxlv = 0.d0
    fluxg = 0.d0
    fluxgh = 0.d0
    fluxgv = 0.d0
    
      
      
       ! We need to calculate the "diffusion" constant D at the interface;
       ! it is defined at the cell centers.  We use the harmonic mean of the
       ! values from the two cells.
       
       !     D1 = perm1 / viscosity_loc_p(jm1)
       !     D2 = perm2 / viscosity_loc_p(jm2)
       !     Dq = (D1 * D2) / (dd2*D1 + dd1*D2)


 
!   calculate average permeability
    D1 = perm1 !* viscosity_loc_p(m2) 
    D2 = perm2 !* viscosity_loc_p(m1) 
    den = dd2*D1 + dd1*D2
    Dq = (perm1 * perm2) / den   ! average k


!****************** for liquid Phase Flux***********************************
    if(((1.D0-SSATG_LOC(m1))>grid%swir(int(icap_loc_p(m1)))).or. &
              ((1.D0-SSATG_LOC(m2))>grid%swir(int(icap_loc_p(m2)))))then
      j=1             ! index to handle varibles defined with context da_nphase
      jm1= j + (m1-1) * grid%nphase
      jm2= j + (m2-1) * grid%nphase
                      

      upweight=0.5d0  ! weight for m1
      if((1.D0-SSATG_LOC(m1)) < grid%swir(int(icap_loc_p(m1))))then
        upweight=0.d0
      else if((1.D0-SSATG_LOC(m2)) < grid%swir(int(icap_loc_p(m2))))then
        upweight=1.d0
      endif
      density_ave = upweight*ddensity_loc_p(jm1)+ &
                     (1.D0-upweight)*ddensity_loc_p(jm2)  

!     gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)
      gravity = (upweight*ddensity_loc_p(jm1)*avgmw_loc_p(jm1) + &
              (1.D0-upweight)*ddensity_loc_p(jm2)*avgmw_loc_p(jm2)) &
              * grid%gravity * grid%delz(nc)

      dphi= -(PPRESSURE_LOC(m2) - PPRESSURE_LOC(m1) & 
                 - pc_loc_p(jm2) + pc_loc_p(jm1)   &
                 - gravity)

      if(dphi>=0.D0)then
        mu=m1
      else
        mu=m2
      end if
      jmu= j + (mu-1) * grid%nphase          

!     print *,'pflow_2ph: ',nc,m1,m2,mu,grid%delz(nc),dphi, &
!     PPRESSURE_LOC(m2),PPRESSURE_LOC(m1),upweight,gravity,ddensity_loc_p(jm1), &
!     ddensity_loc_p(jm2),avgmw_loc_p(jm1),avgmw_loc_p(jm2)

      if((kvr_loc_p(jmu)*Dq)>floweps)then  ! liquid phase is mobile

! calculate mobility with upsteam weight

! Calculate the velocity

!        D1=kvr_loc_p(jm1) 
!        D2=kvr_loc_p(jm2) 
    
!        if(mu==m1)then
!           D0=kvr_loc_p(jm1) 
!        else
!           D0=kvr_loc_p(jm2) 
!        endif
        
         D0=kvr_loc_p(jmu) 
         v_darcy =  Dq * D0 * dphi

 !    store velocities defined at interfaces in PETSc Vec vl at upstream node
      grid%vvl_loc(nc) = v_darcy     ! use for coupling to ptran
      if (n1 > 0) then               ! If the upstream node is not a ghost node...
        vl_p(1+(ip1-1)*grid%nphase+3*grid%nphase*(n1-1)) = v_darcy 
        ! use for print out of velocity
      endif

        q = v_darcy * grid%area(nc)
 !       q = 0.D0 !*********************************************

 !  store velocities defined at interfaces in PETSc Vec vl at upstream node
 !     grid%vvl_loc(nc) = v_darcy     ! use for coupling to ptran
 ! if (n1 > 0) then               ! If the upstream node is not a ghost node...
 !    vl_p(ip1+3*(n1-1)) = v_darcy ! use for print out of velocity
 !   endif
 
       fluxl = fluxl + density_ave * q
       fluxlh = fluxlh +  density_ave * q * hh_loc_p(jmu)
       fluxlv = fluxlv +  density_ave * q *(1.D0-Hen_loc_p(2+(j-1)*grid%nspec+ &
                 (mu-1)*grid%nphase*grid%nspec)*CCONC_LOC(mu))
 
    end if
 end if

!print *,' build r ::convc liq   ',nc,  fluxl,fluxlh ,fluxlv
!************* for Gas phase Flux************************************************* 
if((SSATG_LOC(m1)>eps).or.(SSATG_LOC(m2)>eps))then
    j=2             ! index to handle varibles defined with context da_nphase
    jm1= j + (m1-1) * grid%nphase
    jm2= j + (m2-1) * grid%nphase
                  
    
    upweight=0.5d0  ! weight for m1
    if(SSATG_LOC(m1)<eps)then
       upweight=0.D0
    else if(SSATG_LOC(m2)<eps)then
       upweight=1.D0
    endif
    density_ave = upweight*ddensity_loc_p(jm1)+ &
                  (1.D0-upweight)*ddensity_loc_p(jm2)  
 

    gravity = (upweight*ddensity_loc_p(jm1)*avgmw_loc_p(jm1) + &
              (1.D0-upweight)*ddensity_loc_p(jm2)*avgmw_loc_p(jm2)) &
              * grid%gravity * grid%delz(nc)
 
    dphi= -(PPRESSURE_LOC(m2) - PPRESSURE_LOC(m1)- gravity)

    if(dphi>=0.D0)then
       mu=m1
    else
       mu=m2
    end if
    jmu= j + (mu-1) * grid%nphase        

    if((kvr_loc_p(jmu)*Dq)>floweps)then
      
    !   D1=krwk_loc_p(m1) / viscosity_loc_p(jm1) 
     !  D2=krwk_loc_p(m2) / viscosity_loc_p(jm2)
    
      ! if(mu==m1)then
       !   D0=D1
      ! else
      !    D0=D2
      ! endif

       D0=kvr_loc_p(jmu)
       v_darcy =  Dq * D0 * dphi
! store velocities defined at interfaces in PETSc Vec vl at upstream node
      grid%vvg_loc(nc) = v_darcy     ! use for coupling to ptran
      if (n1 > 0) then               ! If the upstream node is not a ghost node...
        vl_p(2+(ip1-1)*grid%nphase +3*grid%nphase*(n1-1)) = v_darcy 
        ! use for print out of velocity
      endif

       q = v_darcy * grid%area(nc)
!      q = 0.D0 !*****************************************************
       fluxg = fluxg + density_ave * q
       fluxgh = fluxgh + density_ave * q * hh_loc_p(jmu)
       fluxgv = fluxgv + density_ave * q * (1.D0 - CCONC_LOC(mu))
    end if
  !print *,' build r ::convc gas   ',nc,  fluxg,fluxgh ,fluxgv
 end if

 ! heat & tracer residual 
       i1 = ithrm_loc_p(m1)
       i2 = ithrm_loc_p(m2)
       D1 = grid%ckwet(i1)
       D2 = grid%ckwet(i2)

 !heat conductance may involve in saturation later
       Dk = (D1 * D2) / (dd2*D1 + dd1*D2)

       cond = Dk * grid%area(nc)
       hflx = cond * (TTEMP_LOC(m2) - TTEMP_LOC(m1))
       fluxe = fluxlh + fluxgh - hflx

  
! now for simplicity, the difaq is a const and difga is estimated based on it. 
       por1 = porosity_loc_p(m1)
       por2 = porosity_loc_p(m2)
 !      dif1 = df_loc_p(m1)
 !      dif2 = df_loc_p(m2)
       difff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
       diffl=0.D0; diffg=0.D0
      if(((1.D0-SSATG_LOC(m1))>eps).and.((1.D0-SSATG_LOC(m2))>eps))then 
       diffl=difff*0.5d0*(2.d0-SSATG_LOC(m1)-SSATG_LOC(m2))
      endif

      if((SSATG_LOC(m1)>eps).and.(SSATG_LOC(m2))>eps)then 
         difff = (por1 * por2) / (dd2*por1 + dd1*por2) 
        difff = difff * grid%area(nc)*0.25d0
        diffg=(grid%cdiff(int(ithrm_loc_p(m1)))+ &
        grid%cdiff(int(ithrm_loc_p(m2))))*0.5D0
       
        difff = difff *diffg
        diffg=difff*(SSATG_LOC(m1)+SSATG_LOC(m2))* &
        (ddensity_loc_p(2+(m1-1)*grid%nphase)+ ddensity_loc_p(2+(m2-1)*grid%nphase))           
      endif
  !print *,'Res, difff, diffg=',difff,diffg
       difluxg = diffg *  (-CCONC_LOC(m2) + CCONC_LOC(m1))

       j=1
       cw1 =1.D0-(CCONC_LOC(m1))*Hen_loc_p(2+(j-1)*grid%nspec+ &
                  (m1-1)*grid%nphase*grid%nspec)
       cw2 =1.D0-(CCONC_LOC(m2))*Hen_loc_p(2+(j-1)*grid%nspec+ &
                  (m2-1)*grid%nphase*grid%nspec)
     !  print *,'Residual diff cw',cw1,cw2
       difluxl = diffl * grid%area(nc) * (cw2 - cw1)*0.5d0 * &
                  (ddensity_loc_p(1+(m1-1)*grid%nphase)+ &
                 ddensity_loc_p(1+(m2-1)*grid%nphase))

  
  ! difluxl =0.D0; difluxg =0.D0 ! Zero SETTing
   
       fluxcl =  fluxlv - difluxl
       fluxcg =  fluxgv - difluxg
   !fluxcl =0.D0 ; fluxcg =0.D0 ! Zero SETTing
   !    print *,'Residual cf ',fluxg,fluxl,fluxgv,fluxlv, fluxcg,fluxcl
    ! interface on rhs
       if (n1 > 0) then  ! If upstream node is not a ghost node...
          r_p(p1) = r_p(p1) + fluxg + fluxl
          r_p(t1) = r_p(t1) + fluxe            ! heat

          r_p(s1) = r_p(s1) + fluxcl + fluxcg

          if(iphase_loc_p(m1)==6)then
             r_p(c1) = r_p(c1) + fluxcg          ! tracer
          end if
  !        else
  !           r_p(c1) = r_p(c1) + fluxcl
  !        endif
  !     if(iphase_loc_p(m1)==6) r_p(s1) = r_p(s1) + fluxcl + fluxcg
             
       endif

    ! interface on lhs
       if (n2 > 0) then ! If downstream node is not a ghost node...
          r_p(p2) = r_p(p2) - fluxg - fluxl
          r_p(t2) = r_p(t2) - fluxe

          r_p(s2) = r_p(s2) - fluxcl - fluxcg

         if(iphase_loc_p(m2)==6) then
             r_p(c2) = r_p(c2) - fluxcg
          end if
 !        else
 !           r_p(c2) = r_p(c2) - fluxcl
 !         end if
  !      if(iphase_loc_p(m2)==6)  r_p(s2) = r_p(s2) - fluxcl - fluxcg
       endif
  
 ! print *,' build r ::flux c1   ',r_p(c1), r_p(s1),fluxcg,fluxcl,difluxl
 end do
 
 
!*************** Handle boundary conditions*************
!   print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!  print *,'2ph bc-sgbc', grid%myrank, grid%sgbc    
  
  do nc = 1, grid%nconnbc

       m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       ng = grid%nL2G(m)

       if(ng<=0)then
         print *, "Wrong boundary node index... STOP!!!"
         stop
       end if
  
       p1 = 1 + (m-1) * grid%ndof
       t1 = p1 + 1
       c1 = t1 + 1
       s1 = c1 + 1

       ibc = grid%ibconn(nc)
       ip1 = grid%ipermbc(nc)

       select case(ip1)
         case(1)
           perm1 = perm_xx_loc_p(ng)
         case(2)
           perm1 = perm_yy_loc_p(ng)
         case(3)
           perm1 = perm_zz_loc_p(ng)
       end select

       select case(grid%ibndtyp(ibc))
          
          case(2)
          ! solve for pb from Darcy's law given qb /= 0
             grid%pressurebc(:,nc)=PPRESSURE_LOC(ng)
             grid%tempbc(nc) = TTEMP_LOC(ng)
             grid%sgbc(nc) = SSATG_LOC(ng)
             grid%concbc(nc) = CCONC_LOC(ng)
          case(3) 
             grid%tempbc(nc) = TTEMP_LOC(ng)
             grid%sgbc(nc) = SSATG_LOC(ng)
             grid%concbc(nc) = CCONC_LOC(ng)
          end select

! print *,'2ph bc',grid%myrank,nc,m,ng,ibc,grid%ibndtyp(ibc),grid%pressurebc(:,ibc), &
! grid%tempbc(ibc),grid%sgbc(ibc),grid%concbc(ibc),grid%velocitybc(:,ibc)

!   if(grid%ibndtyp(ibc) == 1) then

      !need specify injection phase ratio,conc and pressure
   !   grid%ibndphaseRate(ibc) 
   !   grid%ibndconc(ibc)    ! 
   !   grid%tempbc(ibc)      !1 elements 
   !   grid%pressurebc(ibc)  !nphase elements
!      endif
   
    iicap=  icap_loc_p(ng)     
       
!      print *,'pflow_2pha_bc: ',grid%myrank,' nc= ',nc,' m= ',m, &
!      ' ng= ',ng,' ibc= ',ibc,ip1,iicap, &
!      grid%nconnbc,grid%ibndtyp(ibc),grid%concbc(nc)
     
!   print *,'pflow_2pha-bc: ',ibc,grid%ideriv,grid%ibndtyp(ibc),grid%density_bc,&
!   grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
       
    if (grid%ideriv == 1) then
      call mixture_eos(grid%pressurebc(2,nc),grid%tempbc(nc),&
      grid%concbc(nc),grid%sgbc(nc),&
      grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
      iicap,grid%swir(iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), & !use node's value
      grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,grid%density_bc,&
      grid%d_p_bc,grid%d_t_bc,grid%d_c_bc,grid%d_s_bc,&
      grid%avgmw_bc,grid%avgmw_c_bc,&
      grid%hh_bc,grid%h_p_bc,grid%h_t_bc,grid%h_c_bc,grid%h_s_bc,&
      grid%uu_bc,grid%u_p_bc,grid%u_t_bc,grid%u_c_bc,grid%u_s_bc,&
      grid%df_bc,grid%df_p_bc,grid%df_t_bc,grid%df_s_bc,grid%df_c_bc,&
      grid%hen_bc,grid%hen_p_bc,grid%hen_t_bc,grid%hen_s_bc,grid%hen_c_bc,&
      grid%xxphi_co2_bc(nc),&
      grid%pc_bc,grid%pc_p_bc,grid%pc_t_bc,grid%pc_c_bc,grid%pc_s_bc,&
      grid%kvr_bc,grid%kvr_p_bc,grid%kvr_t_bc,grid%kvr_c_bc,grid%kvr_s_bc,&
      ierr,grid%itable)
    else    
      call mixture_eos_noderiv (grid%pressurebc(2,nc),grid%tempbc(nc),&
      grid%concbc(nc),grid%sgbc(nc),&
      grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
      iicap, grid%swir(iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), &
      grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
      grid%density_bc,grid%avgmw_bc,grid%hh_bc,grid%uu_bc,&
      grid%df_bc,grid%hen_bc,grid%xxphi_co2_bc(nc),grid%pc_bc,grid%kvr_bc,&
      ierr,grid%itable)
    end if

    !   print *, ' boundary index', nc,ng,ibc,grid%ibndtyp(ibc)
    !   print *, ' P  T   C   S  ', grid%pressurebc(1,ibc),grid%tempbc(ibc), &
    !                               grid%concbc(ibc),grid%sgbc(ibc)
    !   print *,' hh,den   ',grid%hh_bc(1:2),grid%density_bc(1:2)

!print *,' Gotten BC properties ', ibc,grid%ibndtyp(ibc),iicap
!print *,grid%pressurebc(2,ibc),grid%tempbc(ibc),grid%concbc(ibc),grid%sgbc(ibc)
!print *,grid%density_bc,grid%avgmw_bc
!print *,grid%hh_bc,grid%uu_bc,grid%df_bc,grid%hen_bc,grid%pc_bc,grid%kvr_bc

   flux = 0.d0
   fluxh = 0.d0
   fluxv = 0.d0
   
   fluxl = 0.d0
   fluxlh = 0.d0
   fluxlv = 0.d0

   fluxg = 0.d0
   fluxgh = 0.d0
   fluxgv = 0.d0

   select case (grid%ibndtyp(ibc))
     case(1)  ! Dirichlet BC for p, T, C, S...
       !***********Liquid phase ************************
       Dq = perm1 / grid%distbc(nc)
       if(((1.D0-grid%sgbc(nc)) > grid%swir(int(icap_loc_p(ng)))).or.&
           ((1.D0-SSATG_LOC(ng)) > grid%swir(int(icap_loc_p(ng))))) then
          j=1
         ! index to handle varibles defined with context da_nphase
          jng = j + (ng-1) * grid%nphase             
          
          upweight=1.d0  ! weight for virtual BC node
          if((1.D0-grid%sgbc(nc)) < grid%swir(int(icap_loc_p(ng))))then
            upweight=0.d0
          endif
          density_ave = upweight* grid%density_bc(j)+ &
                        (1.D0-upweight)*ddensity_loc_p(jng)  


! if liquid phase injected, use the krw in the boundary node
          
          gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
              (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
              * grid%gravity * grid%delzbc(nc)
          dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc) & 
                - pc_loc_p(jng) + grid%pc_bc(j)   &
                - gravity)
        
          if(dphi >= 0.D0)then
             !upstream is BC
             ukvr=grid%kvr_bc(j)
             uhh=grid%hh_bc(j)
             uconc=1.D0-grid%concbc(nc)*grid%hen_bc(2)
          else
             !upstream is node
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
             uconc=1.D0-CCONC_LOC(ng)*hen_loc_p(i2_hencoeff)
          end if
      

          if((ukvr*Dq) > floweps)then
! note: darcy vel. is positive for flow INTO boundary node
            v_darcy =  Dq *  ukvr  * dphi
            grid%vvlbc(nc) = v_darcy    
            if (m > 0) then    ! If the upstream node is not a ghost node...
              vl_p(1+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
              ! use for print out of velocity
            endif

            q = v_darcy * grid%areabc(nc)
            !q = 0.d0
! the original THCResidual is not consistant here in the averaging aspect
! So.. changed
        
            fluxl = fluxl + q * density_ave
            fluxlh = fluxlh +  density_ave * q * uhh
            fluxlv = fluxlv +  density_ave * q * uconc
!           print *,'BC1 liq f ',dphi,PPRESSURE_LOC(ng),grid%pressurebc(2,ibc),&
!                  pc_loc_p(jng),grid%pc_bc(j)

 !          print *,'          ', grid%sgbc(ibc),SSATG_LOC(ng),&
 !          PPRESSURE_LOC(ng)- pc_loc_p(jng),grid%pressurebc(2,ibc)-grid%pc_bc(j)
 !          print *,'          ',fluxl, fluxlh ,fluxlv
          end if
!         print *,'pflow_2ph_bc1: ',nc,q,uconc,density_ave,floweps
       end if

!************** Gas Phase injection************************
          if((grid%sgbc(nc)>eps).or.(SSATG_LOC(ng)>eps))then
            j=2     ! index to handle varibles defined with context da_nphase
            jng= j + (ng-1) * grid%nphase
 
            upweight=1.d0  ! weight for virtual BC node
            if(grid%sgbc(nc)<eps)then
              upweight=0.d0
            endif
            density_ave = upweight* grid%density_bc(j)+ &
                          (1.D0-upweight)*ddensity_loc_p(jng)  

            gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                       (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                       * grid%gravity * grid%delzbc(nc)
            dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc)- gravity)

              
          if(dphi>=0.D0)then
             !upstream is BC
             ukvr=grid%kvr_bc(j)
             uhh=grid%hh_bc(j)
             uconc=1.D0-grid%concbc(nc)
          else
             !upstream is node
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc=1.D0-CCONC_LOC(ng)
            
          end if
          if((ukvr*Dq)>floweps)then
            ! note: darcy vel. is positive for flow INTO boundary node
            v_darcy =  Dq *  ukvr  * dphi
            grid%vvgbc(nc) = v_darcy
            if (m > 0) then      ! If the upstream node is not a ghost node...
              vl_p(2+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
              ! use for print out of velocity
            endif
    
            q = v_darcy * grid%areabc(nc)
            !q = 0.d0
! the original THCResidual is not consistant here in the averaging aspect
! So.. changed
        
            fluxg = fluxg+ q * density_ave
            fluxgh = fluxgh +  density_ave * q * uhh
            fluxgv = fluxgv +  density_ave * q * uconc
  ! print *,'BC1 gas f ',dphi,PPRESSURE_LOC(ng),grid%pressurebc(2,ibc),fluxg, fluxgv
                 

   !          print *,'          ', grid%sgbc(ibc),SSATG_LOC(ng)
               
          end if
       end if
 
      r_p(p1) = r_p(p1) - fluxg - fluxl
             
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      r_p(t1)=r_p(t1) - fluxgh - fluxlh + cond * &
      (TTEMP_LOC(ng) - grid%tempbc(nc))
                     

          
! Diffusion term           
      diffl=0.D0; diffg=0.D0
      difff =grid%difaq * porosity_loc_p(ng)/ grid%distbc(nc)
      if(((1.D0-grid%sgbc(nc))>eps).and.((1.D0-SSATG_LOC(ng))>eps))then 
        diffl=difff*0.5d0*(2.d0-grid%sgbc(nc)-SSATG_LOC(ng))
      endif
      if((grid%sgbc(nc)>eps).and.(SSATG_LOC(ng))>eps)then 
        diffg=grid%cdiff(int(ithrm_loc_p(ng)))*porosity_loc_p(ng)/ grid%distbc(nc)
        diffg=diffg*0.5d0*(grid%sgbc(nc)+SSATG_LOC(ng))
      endif

      difluxg = diffg * grid%areabc(nc) * (-CCONC_LOC(ng)+ &
                  grid%concbc(nc))*0.5d0*( grid%density_bc(2)+&
                  ddensity_loc_p(2+(ng-1)*grid%nphase))

      i1_hencoeff=2
      cw1 =1.D0-(grid%concbc(nc))*grid%hen_bc(i1_hencoeff)
      i2_hencoeff=2 + (ng-1)*grid%nphase*grid%nspec
      cw2 =1.D0-(CCONC_LOC(ng))*hen_loc_p(i2_hencoeff)

      difluxl = diffl * grid%areabc(nc) * 0.5d0 * &
      (grid%density_bc(1)+ ddensity_loc_p(1+(ng-1)*grid%nphase))*(cw2 - cw1)
             
          
      fluxcl =  fluxlv - difluxl
      fluxcg =  fluxgv - difluxg
          
      r_p(s1)=r_p(s1) - fluxcl - fluxcg
             
      if(iphase_loc_p(ng)==6)then
        r_p(c1)=r_p(c1) - fluxcg
      endif         
       
!             if(iphase_loc_p(ng)==6)   r_p(s1)=r_p(s1) - fluxcl - fluxcg
! print *,'BC1 c,diffl  ', grid%concbc(ibc),difluxl,difluxg,fluxcl,fluxcg
! print *,'             ',  fluxl,fluxg,grid%sgbc(ibc),SSATG_LOC(ng)                                    

!---------------------------------------------------------------------------

    case(2)    ! Constant velocity q, grad T, C = 0
 
! Notice :: different from THCResidual, that if q>0, should specify injection 
! here we control the total flux: include contribution from diffusion if q>0
! so there is no diffuaion term appeared
 if((dabs(grid%velocitybc(1,nc))+dabs(grid%velocitybc(2,nc)))>floweps)then
       do j=1,grid%nphase
          jng = j + (ng-1) * grid%nphase
          fluxv=0.D0; flux=0.D0;  fluxh=0.D0
          v_darcy = grid%velocitybc(j,nc)
      
      select case(j)
      case(1)
        grid%vvlbc(nc) = v_darcy
            if (m > 0) then          ! If the upstream node is not a ghost node...
              vl_p(1+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
              ! use for print out of velocity
            endif
  
      case(2)
        grid%vvgbc(nc) = v_darcy 
            if (m > 0) then          ! If the upstream node is not a ghost node...
              vl_p(2+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
              ! use for print out of velocity
            endif

           end select  
          
      if(v_darcy >0.d0)then 
             q = v_darcy * grid%density_bc(j) * grid%areabc(nc)
             !q = 0.d0
             flux = flux - q
             fluxh = fluxh - q  * grid%hh_bc(j) 
             select case(j)
             case(1) 
                i1_hencoeff=2+(j-1)*grid%nspec
                cw =1.D0-(grid%concbc(nc))*grid%hen_bc(i1_hencoeff)
                fluxv = fluxv - q  * cw
             case(2) 
                fluxv = fluxv - q  * (1.D0-grid%concbc(nc))
              end select
! add diffusion and heat conduction term here
          else
             q = v_darcy * ddensity_loc_p(jng) * grid%areabc(nc)
             !q = 0.d0
             flux = flux - q
             fluxh = fluxh -q * hh_loc_p(jng)
             select case(j)
             case(2) 
                fluxv = fluxv - q * (1.D0-CCONC_LOC(ng))
             case(1)  
                i2_hencoeff=2+(ng-1)*grid%nphase*grid%nspec
                cw =1.D0-(CCONC_LOC(ng))* hen_loc_p(i2_hencoeff)
                 fluxv = fluxv - q  * cw
            end select
          end if
       if((j==2).and.(SSATG_LOC(ng)>eps)) r_p(c1)=r_p(c1)-fluxv 
       if((j==1).and.(SSATG_LOC(ng)<=eps)) r_p(c1)=r_p(c1)-fluxv 
!
       r_p(p1) = r_p(p1) - flux
       r_p(t1) = r_p(t1) - fluxh
        if(iphase_loc_p(m2)==6)  r_p(s1) = r_p(s1) - fluxv
    end do
 endif
 
  case(3) ! fixed p, grad T, C = 0
         Dq = perm1 / grid%distbc(nc)
      
      if((1.D0-SSATG_LOC(ng))>grid%swir(int(icap_loc_p(ng))))then
          j=1
         ! index to handle varibles defined with context da_nphase
          jng = j + (ng-1) * grid%nphase             
          
          upweight=0.D0  ! weight for virtual BC node
!          if((1.D0-grid%sgbc(ibc))<grid%swir(int(icap_loc_p(ng))))then
!             upweight=0.
!          endif
          density_ave = upweight* grid%density_bc(j)+ &
                        (1.D0-upweight)*ddensity_loc_p(jng)  


! if liquid phase injected, use the krw in the boundary node
          
          gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
              (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
              * grid%gravity * grid%delzbc(nc)
          dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc) & 
                 - pc_loc_p(jng) + grid%pc_bc(j)   &
                 - gravity)

     
     
     !     dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(1,ibc) - gravity)
        
     !     if(dphi>0)then
             !upstream is BC
     !        ukvr=grid%kvr_bc(j)
     !        uhh=grid%hh_bc(j)
     !        i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
     !        uconc=1.D0-(1.D0-CCONC_LOC(ng))*hen_loc_p(i2_hencoeff)
     !       else
             !upstream is node
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
             uconc=1.D0-(CCONC_LOC(ng))*hen_loc_p(i2_hencoeff)
    !      end if
      

          if((ukvr*Dq)>floweps)then
! note: darcy vel. is positive for flow INTO boundary node
            v_darcy =  Dq *  ukvr  * dphi
            grid%vvlbc(nc) = v_darcy
            if (m > 0) then   ! If the upstream node is not a ghost node...
              vl_p(1+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
              ! use for print out of velocity
            endif
              
             q = v_darcy * grid%areabc(nc)
             !q = 0.d0
! the original THCResidual is not consistent here in the averaging aspect
! So.. changed
        
             fluxl = fluxl + q * density_ave
             fluxlh = fluxlh +  density_ave * q * uhh
             fluxlv = fluxlv +  density_ave * q * uconc
!             print *,'BC 3',fluxl,(CCONC_LOC(ng)) ,uconc,hen_loc_p(i2_hencoeff)
          end if
       end if
!************** Gas Phase injection************************
          if(SSATG_LOC(ng)>eps)then
             j=2        ! index to handle varibles defined with context da_nphase
             jng= j + (ng-1) * grid%nphase
 
             upweight=0.D0  ! weight for virtual BC node
        
             density_ave = upweight* grid%density_bc(j)+ &
                          (1.D0-upweight)*ddensity_loc_p(jng)  

             gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                       (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                       * grid%gravity * grid%delzbc(nc)
             dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc)- gravity)

              
     !     if(dphi>0)then
             !upstream is BC
      !       ukvr=grid%kvr_bc(j)
      !       uhh=grid%hh_bc(j)
      !       uconc=grid%concbc(ibc)
      !    else
             !upstream is node
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc=1.D0-CCONC_LOC(ng)
            
       !   end if
        
             if((ukvr*Dq)>floweps)then
               ! note: darcy vel. is positive for flow INTO boundary node
               v_darcy =  Dq *  ukvr  * dphi
         grid%vvgbc(nc) = v_darcy
               if (m > 0) then   ! If the upstream node is not a ghost node...
                 vl_p(2+(ip1-1)*grid%nphase +3*grid%nphase*(m-1)) = -v_darcy 
                 ! use for print out of velocity
               endif
              
               q = v_darcy * grid%areabc(nc)
               !q = 0.d0
! the original THCResidual is not consistant here in the averaging aspect
! So.. changed
        
               fluxg = fluxg+ q * density_ave
               fluxgh = fluxgh +  density_ave * q * uhh
               fluxgv = fluxgv +  density_ave * q * uconc
             end if
           end if
     
            r_p(p1) = r_p(p1) - fluxg - fluxl
            r_p(t1)=r_p(t1) - fluxgh - fluxlh

            r_p(s1)=r_p(s1) - grid%fc * fluxlv -grid%fc * fluxgv

            if(iphase_loc_p(m2)==6)then
               r_p(c1)=r_p(c1) - grid%fc * fluxgv 
            end if

         end select
     
!      print *,'pflow_2pha-bc: ',nc,ibc,grid%ibndtyp(ibc), &
!      grid%density_bc(1),grid%density_bc(2),&
!      grid%pressurebc(2,nc),grid%tempbc(nc),grid%concbc(nc),grid%sgbc(nc)
  
      end do

  do n = 1, grid%nlmax  ! For each local node do switch equation order 3<-->4
    ng = grid%nL2G(n) 
    p1 = 1 + (n-1)*grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1
    s1 = c1 + grid%npricomp

   ! if(SSATG_LOC(ng)<=eps)then
    uhh= r_p(c1)
    r_p(c1)=r_p(s1)
    r_p(s1)=uhh
  ! end if
  end do
 call VecRestoreArrayF90(grid%accum, accum_p, ierr) 
 call VecRestoreArrayF90(r, r_p, ierr)
 call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
 call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
 call VecRestoreArrayF90(grid%avgmw_loc, avgmw_loc_p, ierr)
!  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(grid%uu, uu_p, ierr)
! call VecRestoreArrayF90(grid%h, h_p, ierr)
 
  call VecRestoreArrayF90(grid%df_loc, df_loc_p, ierr)
  call VecRestoreArrayF90(grid%hen_loc, hen_loc_p, ierr)
  call VecRestoreArrayF90(grid%pcw_loc, pc_loc_p, ierr)
  call VecRestoreArrayF90(grid%kvr_loc, kvr_loc_p, ierr)

  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
 
  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(grid%vl, vl_p, ierr)
  call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)
 

  if (grid%ideriv == 1) then
   call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_c_loc, d_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_s_loc, d_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%avgmw_c,avgmw_c_p, ierr)
    call VecRestoreArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_c_loc, h_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_s_loc, h_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%u_p, u_p_p, ierr)
    call VecRestoreArrayF90(grid%u_t, u_t_p, ierr)
    call VecRestoreArrayF90(grid%u_c, u_c_p, ierr)
    call VecRestoreArrayF90(grid%u_s, u_s_p, ierr)
    call VecRestoreArrayF90(grid%df_p_loc, df_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_t_loc, df_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_c_loc, df_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_s_loc, df_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_p_loc, hen_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_t_loc, hen_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_c_loc, hen_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_s_loc, hen_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_p_loc, pc_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_t_loc, pc_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_c_loc, pc_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_s_loc, pc_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_p_loc, kvr_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_t_loc, kvr_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_c_loc, kvr_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_s_loc, kvr_s_loc_p, ierr)
  endif

!print *,'xxxxxxxxx ::...........'; call VecView(xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
!print *,'Residual ::...........'; call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

 end subroutine TTPHASEResidual
 

! --------------------------------------------------------------------- 


  subroutine TTPHASEJacobin(snes,xx,A,B,flag,grid,ierr)
      use mixture_module  
      use water_eos_module
  
    implicit none

    SNES, intent(in) :: snes
    Vec, intent(in) :: xx
    Mat, intent(inout) :: A, B
    type(pflowGrid), intent(inout) :: grid
   ! integer, intent(inout) :: flag
    MatStructure flag

    integer :: ierr
    integer :: n, ng, nc
    integer :: i1, i2, j, jn, jng, jm1, jm2,jmu
    integer :: m, m1, m2, mu, n1, n2, ip1, ip2 
    integer :: i1_hencoeff, i2_hencoeff,i1_hencoeff_dc,i2_hencoeff_dc, &
               iu_hencoeff, iu_hencoeff_dc
    integer :: p1,p2,t1,t2,c1,c2,s1,s2
    integer :: ibc  ! Index that specifies a boundary condition block.
    real*8 ::  v_darcy, q

    PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), &
               ddensity_loc_p(:),&
               phis_p(:),  &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
               
!   PetscScalar, pointer :: d_p_p(:), d_p_loc_p(:), &
!              d_t_p(:), d_t_loc_p(:), &
!              d_c_p(:), d_c_loc_p(:), &
!              d_s_p(:), d_s_loc_p(:)
               
    PetscScalar, pointer :: d_p_loc_p(:), &
               d_t_loc_p(:), &
               d_c_loc_p(:), &
               d_s_loc_p(:)
               
               
!   PetscScalar, pointer :: avgmw_p(:),avgmw_loc_p(:),avgmw_c_p(:),avgmw_c_loc_p(:),&
!              h_p(:),  hh_p(:), hh_loc_p(:), &
!              h_p_p(:), h_p_loc_p(:), &
!              h_t_p(:), h_t_loc_p(:), &
!              h_c_p(:), h_c_loc_p(:), &
!              h_s_p(:), h_s_loc_p(:)
               
    PetscScalar, pointer :: avgmw_loc_p(:),avgmw_c_loc_p(:),&
               hh_loc_p(:), &
               h_p_loc_p(:), &
               h_t_loc_p(:), &
               h_c_loc_p(:), &
               h_s_loc_p(:)
               
               
    PetscScalar, pointer :: uu_p(:),  &
               u_p_p(:),u_t_p(:),u_c_p(:),u_s_p(:)
               
               
!   PetscScalar, pointer :: hen_p(:),  hen_loc_p(:), &  
!              hen_p_p(:), hen_p_loc_p(:), &
!              hen_t_p(:), hen_t_loc_p(:), &
!              hen_c_p(:), hen_c_loc_p(:), &
!              hen_s_p(:), hen_s_loc_p(:)
               
    PetscScalar, pointer :: hen_loc_p(:), &  
               hen_p_loc_p(:), &
               hen_t_loc_p(:), &
               hen_c_loc_p(:), &
               hen_s_loc_p(:)
               
               
!   PetscScalar, pointer :: df_p(:), df_loc_p(:), &
!              df_p_p(:), df_p_loc_p(:),&
!              df_t_p(:), df_t_loc_p(:),&
!              df_c_p(:), df_c_loc_p(:),&
!              df_s_p(:), df_s_loc_p(:)
               
    PetscScalar, pointer :: df_loc_p(:), &
               df_p_loc_p(:),&
               df_t_loc_p(:),&
               df_c_loc_p(:),&
               df_s_loc_p(:)
               
! PetscScalar, pointer :: pc_p(:), pc_loc_p(:),&
!                         pc_p_p(:), pc_p_loc_p(:),&
!                         pc_t_p(:), pc_t_loc_p(:),&
!                         pc_c_p(:), pc_c_loc_p(:),&
!                         pc_s_p(:), pc_s_loc_p(:),&
!                         kvr_p(:), kvr_loc_p(:),&
!                         kvr_p_p(:), kvr_p_loc_p(:),&
!                         kvr_t_p(:), kvr_t_loc_p(:),&
!                         kvr_c_p(:), kvr_c_loc_p(:),&
!                         kvr_s_p(:), kvr_s_loc_p(:)
               
  PetscScalar, pointer :: pc_loc_p(:),&
                          pc_p_loc_p(:),&
                          pc_t_loc_p(:),&
                          pc_c_loc_p(:),&
                          pc_s_loc_p(:),&
                          kvr_loc_p(:),&
                          kvr_p_loc_p(:),&
                          kvr_t_loc_p(:),&
                          kvr_c_loc_p(:),&
                          kvr_s_loc_p(:)


  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  integer :: iicap,ii,jj,iiphas,iiphas1,iiphas2
  integer ibc_hencoeff
  real*8 :: cond, gravity, SSATW, acc, &
            density_ave, voldt, pvoldt
  real*8 :: fluxl, fluxlh, fluxlv, fluxg, fluxgh, fluxgv, &
            flux, fluxh, fluxv, difff, diffg, diffl
  real*8 :: dd1, dd2, dd, f1, f2, den
! real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2
  real*8 :: por1, por2, perm1, perm2
  real*8 :: qu_rate, p_vapor,sat_pressure_t
! real*8 :: cg1,cg2,cg,cg_p,cg_t,cg_s,cg_c
  real*8 :: Dk, Dq,D0, Dphi, gdz  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: xxlw,xxla,xxgw,xxga,cw,cw1,cw2,cwu, sat_ave
  real*8 :: ra(1:5,1:8),tempvar(1:8), devq(4,2)  
  real*8 :: uhh, uconc, ukvr
  real*8 :: upweight,m1weight,m2weight,mbweight,mnweight
  real*8 :: cw1_p,cw1_t,cw1_c,cw1_s,cw2_p,cw2_t,cw2_c,cw2_s,cw_p,cw_t,cw_c,cw_s
  real*8 :: blkmat11(1:4,1:4),  blkmat12(1:4,1:4), blkmat21(1:4,1:4), blkmat22(1:4,1:4)

!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
 flag = SAME_NONZERO_PATTERN

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
 
 call DAGlobalToLocalBegin(grid%da_1_dof, grid%iphas, &
                           INSERT_VALUES, grid%iphas_loc, ierr)
 call DAGlobalToLocalEnd(grid%da_1_dof, grid%iphas, &
                          INSERT_VALUES, grid%iphas_loc, ierr)

  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)

! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%avgmw_loc, avgmw_loc_p, ierr)
!  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(grid%uu, uu_p, ierr)
  call VecGetArrayF90(grid%df_loc, df_loc_p, ierr)
  call VecGetArrayF90(grid%hen_loc, hen_loc_p, ierr)
  call VecGetArrayF90(grid%pcw_loc, pc_loc_p, ierr)
  call VecGetArrayF90(grid%kvr_loc, kvr_loc_p, ierr)
 
 !print *,' In 2ph Jacobian :: non got pointers'

    call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecGetArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecGetArrayF90(grid%d_c_loc, d_c_loc_p, ierr)
    call VecGetArrayF90(grid%d_s_loc, d_s_loc_p, ierr)
    call VecGetArrayF90(grid%avgmw_c_loc,avgmw_c_loc_p, ierr)
    call VecGetArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecGetArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecGetArrayF90(grid%h_c_loc, h_c_loc_p, ierr)
    call VecGetArrayF90(grid%h_s_loc, h_s_loc_p, ierr)
    call VecGetArrayF90(grid%u_p, u_p_p, ierr)
    call VecGetArrayF90(grid%u_t, u_t_p, ierr)
    call VecGetArrayF90(grid%u_c, u_c_p, ierr)
    call VecGetArrayF90(grid%u_s, u_s_p, ierr)
    call VecGetArrayF90(grid%df_p_loc, df_p_loc_p, ierr)
    call VecGetArrayF90(grid%df_t_loc, df_t_loc_p, ierr)
    call VecGetArrayF90(grid%df_c_loc, df_c_loc_p, ierr)
    call VecGetArrayF90(grid%df_s_loc, df_s_loc_p, ierr)
    call VecGetArrayF90(grid%hen_p_loc, hen_p_loc_p, ierr)
    call VecGetArrayF90(grid%hen_t_loc, hen_t_loc_p, ierr)
    call VecGetArrayF90(grid%hen_c_loc, hen_c_loc_p, ierr)
    call VecGetArrayF90(grid%hen_s_loc, hen_s_loc_p, ierr)
    call VecGetArrayF90(grid%pc_p_loc, pc_p_loc_p, ierr)
    call VecGetArrayF90(grid%pc_t_loc, pc_t_loc_p, ierr)
    call VecGetArrayF90(grid%pc_c_loc, pc_c_loc_p, ierr)
    call VecGetArrayF90(grid%pc_s_loc, pc_s_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_p_loc, kvr_p_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_t_loc, kvr_t_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_c_loc, kvr_c_loc_p, ierr)
    call VecGetArrayF90(grid%kvr_s_loc, kvr_s_loc_p, ierr)
 !print *,' In 2ph Jacobian ::  got pointers 1'

  
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(grid%iphas_loc, iphase_loc_p, ierr)

  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
  endif
!print *,' In 2ph Jacobian ::  got pointers 2'
! ********************************************************************

! Accumulation terms

  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)   !get ghosted index
    ra=0.D0
   ! Remember, the matrix index starts from (0,0)
    p1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
    t1 = 1 + p1           ! = 2 + (ng-1)*grid%ndof-1
    c1 = 1 + t1           ! = 3 + (ng-1)*grid%ndof-1
    s1=  c1+ grid%npricomp

    j=1   
    jn = j + (n-1)*grid%nphase
    jng = j + (ng-1)*grid%nphase

    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    SSATW=1.D0-SSATG_LOC(ng)
    iiphas=iphase_loc_p(ng)
 ! pressure equation    
    select case (iiphas)
    case(6)
          ra(1,1)=pvoldt*(d_p_loc_p(jng+1)*SSATG_LOC(ng)+d_p_loc_p(jng)*SSATW)
          ra(1,2)=pvoldt *(SSATW * d_t_loc_p(jng)+&
               SSATG_LOC(ng)*d_t_loc_p(jng+1))
          ra(1,3)=pvoldt *(SSATW * d_c_loc_p(jng)+SSATG_LOC(ng)*d_c_loc_p(jng+1))
          ra(1,4)=pvoldt *(-ddensity_loc_p(jng)+ddensity_loc_p(jng+1)+ &
              SSATW*d_s_loc_p(jng) + SSATG_LOC(ng)*d_s_loc_p(jng+1))
    case(2)
            ra(1,1)=pvoldt * d_p_loc_p(jng)
            ra(1,2)=pvoldt * d_t_loc_p(jng)
            ra(1,3)=pvoldt * d_c_loc_p(jng)
            ra(1,4)=0.D0
    case(4)   
            ra(1,1)=pvoldt * d_p_loc_p(jng+1)
            ra(1,2)=pvoldt * d_t_loc_p(jng+1)
            ra(1,3)=pvoldt * d_c_loc_p(jng+1)
            ra(1,4)= 0.D0
    case default    
            print *,'Wrong phase type'
            stop
    end select
 !energy equation
 !   ul=hh_loc_p(jng)- PPRESSURE_LOC(ng)/ volume_p(n)
 !   ug= hh_loc_p(jng+1)- PPRESSURE_LOC(ng)/ volume_p(n)   
   select case (iiphas)
    case(6)
       ra(2,1) = pvoldt * (SSATW*ddensity_loc_p(jng)*u_p_p(jn) &
              + SSATG_LOC(ng)*d_p_loc_p(jng+1)*uu_p(jn+1) &
              + SSATW * d_p_loc_p(jng)* uu_p(jn) &
              + SSATG_LOC(ng)*ddensity_loc_p(jng+1)*u_p_p(jn+1)) 
      !        print *, ' T Accum Jac: ', n, pvoldt, SSATW, ddensity_loc_p(jng),u_p_p(jn), &
      !                    SSATG_LOC(ng),d_p_loc_p(jng+1) ,uu_p(jn+1), &
      !                     d_p_loc_p(jng),  uu_p(jn), &
      !                     ddensity_loc_p(jng+1),u_p_p(jn+1) 
       ! omit dPc/dPg 
       ra(2,2) = pvoldt * (SSATW*d_t_loc_p(jng)*uu_p(jn) + SSATW*&
              ddensity_loc_p(jng)*u_t_p(jn) + SSATG_LOC(ng) *&
              d_t_loc_p(1+jng)*uu_p(jn+1) + SSATG_LOC(ng) *  &
              ddensity_loc_p(jng+1)*u_t_p(jn+1)) + &
              (1.d0-porosity_loc_p(ng))*grid%dencpr(int(ithrm_loc_p(ng))) * &
               voldt  
               ! Rock den and cp  taken as constant
       ra(2,3)=  pvoldt * (SSATW*d_c_loc_p(jng)*uu_p(jn) + SSATW*&
              ddensity_loc_p(jng)*u_c_p(jn) + SSATG_LOC(ng) *&
              d_c_loc_p(1+jng)*uu_p(jn+1) + SSATG_LOC(ng) *  &
              ddensity_loc_p(jng+1)*u_c_p(jn+1))  

       ra(2,4)=  pvoldt * ((SSATW*d_s_loc_p(jng)-ddensity_loc_p(jng))* &
              uu_p(jn)+SSATW*ddensity_loc_p(jng)*u_s_p(jn)+&
              ddensity_loc_p(jng+1)*uu_p(jn+1)+ ddensity_loc_p(jng+1)*&
              SSATG_LOC(ng)*u_s_p(jn+1)+SSATG_LOC(ng)*d_s_loc_p(jng+1)* uu_p(jn+1))
    case(2)
       ra(2,1) = pvoldt * (ddensity_loc_p(jng)*u_p_p(jn) &
            +  d_p_loc_p(jng)* uu_p(jn))
       ra(2,2) = pvoldt * (d_t_loc_p(jng)*uu_p(jn) + &
              ddensity_loc_p(jng)*u_t_p(jn)) + &
              (1.d0-porosity_loc_p(ng))*grid%dencpr(int(ithrm_loc_p(ng))) * &
               voldt  
               ! Rock den and cp  taken as constant
       ra(2,3)=  pvoldt * (d_c_loc_p(jng)*uu_p(jn) + &
              ddensity_loc_p(jng)*u_c_p(jn))  

       ra(2,4)= 0.D0
    case(4)
       ra(2,1) = pvoldt * ( d_p_loc_p(jng+1)*uu_p(jn+1) &
            + ddensity_loc_p(jng+1)*u_p_p(jn+1)) 
       ! omit dPc/dPg 
       ra(2,2) = pvoldt * ( d_t_loc_p(1+jng)*uu_p(jn+1) +   &
              ddensity_loc_p(jng+1)*u_t_p(jn+1)) + &
              (1.d0-porosity_loc_p(ng))*grid%dencpr(int(ithrm_loc_p(ng))) * &
               voldt  
               ! Rock den and cp  taken as constant
       ra(2,3)=  pvoldt * (  d_c_loc_p(1+jng)*uu_p(jn+1)&
            +  ddensity_loc_p(jng+1)*u_c_p(jn+1))  
       ra(2,4)= 0.D0
    end select



    i1_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
    i1_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(ng-1)*grid%nphase*grid%nspec*&
                    grid%npricomp ! +(jc-1)*grid%npricomp
    xxga = CCONC_LOC(ng)
    xxgw = 1.D0-xxga
    xxla = xxga *Hen_loc_p(i1_hencoeff)
    xxlw = 1.D0 - xxla
    cw1_p= - Hen_p_loc_p(i1_hencoeff) * xxga
    cw1_t= - Hen_t_loc_p(i1_hencoeff) * xxga
    cw1_c= - Hen_c_loc_p(i1_hencoeff_dc) * xxga - Hen_loc_p(i1_hencoeff)
    cw1_s= - Hen_s_loc_p(i1_hencoeff) * xxga


 select case (iiphas)
    case(6)
       ra(4,1) = pvoldt * ( SSATW * d_p_loc_p(jng) * xxlw+ &
            SSATG_LOC(ng)* d_p_loc_p(jng+1)* xxgw +   &
            SSATW * ddensity_loc_p(jng) * cw1_p)
       ra(4,2) = pvoldt * ( SSATW * d_t_loc_p(jng) * xxlw+ &
            SSATG_LOC(ng)* d_t_loc_p(jng+1)* xxgw +   &
            SSATW * ddensity_loc_p(jng) *cw1_t )
       ra(4,3) = pvoldt * ( SSATW * d_c_loc_p(jng) * xxlw -&
            SSATG_LOC(ng) * ddensity_loc_p(jng+1)+ &
            SSATG_LOC(ng)* d_c_loc_p(jng+1)* xxgw +  &
            SSATW * ddensity_loc_p(jng) *cw1_c)
       ra(4,4)=  pvoldt * ( -ddensity_loc_p(jng)*xxlw + &
            ddensity_loc_p(jng+1)*xxgw + &
            d_s_loc_p(jng) * SSATW *xxlw + &
            ddensity_loc_p(jng) * SSATW * cw1_s) 
    case(2)
       ra(4,1) = pvoldt * ( d_p_loc_p(jng) * xxlw+ &
             ddensity_loc_p(jng) * cw1_p)
       ra(4,2) = pvoldt * (  d_t_loc_p(jng) * xxlw+ &
             ddensity_loc_p(jng) *cw1_t )
       ra(4,3) = pvoldt * (  d_c_loc_p(jng) * xxlw +&
             ddensity_loc_p(jng) *cw1_c)
       ra(4,4)= 0.D0

    case(4)
       ra(4,1) = pvoldt * d_p_loc_p(jng+1)* xxgw 
       ra(4,2) = pvoldt * d_t_loc_p(jng+1)* xxgw 
       ra(4,3) = pvoldt * ( -  ddensity_loc_p(jng+1)+ d_c_loc_p(jng+1)* xxgw )
       ra(4,4)= 0.D0
     
    end select

!  print *,'energy der dr/dc ',n,ra(2,3),SSATW*d_c_loc_p(jng)*uu_p(jn), &
!            SSATW*ddensity_loc_p(jng)*u_c_p(jn), SSATG_LOC(ng)* d_c_loc_p(1+jng)*uu_p(jn+1),& 
!             SSATG_LOC(ng) *  ddensity_loc_p(jng+1)*u_c_p(jn+1)


 ! concentration equation
   select case (iiphas)
    case(6)
       ra(3,1) = grid%ret * pvoldt * SSATG_LOC(ng) *d_p_loc_p(jng+1)*(1.D0-CCONC_LOC(ng))
       ra(3,2) = grid%ret * pvoldt * SSATG_LOC(ng) *d_t_loc_p(jng+1)*(1.D0-CCONC_LOC(ng))
       ra(3,3) = grid%ret * pvoldt * (d_c_loc_p(jng+1)* (1.D0-CCONC_LOC(ng)) *&
            SSATG_LOC(ng) - ddensity_loc_p(jng+1)*SSATG_LOC(ng))
       ra(3,4) = grid%ret * pvoldt * (ddensity_loc_p(jng+1)*(1.D0-CCONC_LOC(ng))+ &
            SSATG_LOC(ng) *d_s_loc_p(jng+1)*(1.D0-CCONC_LOC(ng)))
    case(2)
        ra(3,1)=0.D0 ; ra(3,2)=0.D0
        ra(3,3)=0.D0 ; ra(3,4)=1.D0
    case(4)
        ra(3,1)=0.D0 ; ra(3,2)=0.D0
        ra(3,3)=0.D0 ; ra(3,4)=1.D0
    end select


!  if(n==1)then
!     print *,'node 1 Jacobian::'
!     print *,' density',ddensity_loc_p(jng:jng+1)
!     print *,'    Xmol Sw', n,xxgw,xxga,xxlw,xxla, SSATW
!      print *, ra(3,4),ra(4,4),ra(5,4)
!   end if
 ! reaction term in concentration equation

 if(iiphas==6)then
  tempvar=0. 
  qu_rate=0.0  ! evaporation get positive
    p_vapor=PPRESSURE_LOC(ng)*xxgw
   
    call psat1(TTEMP_LOC(ng),sat_pressure,sat_pressure_t,ierr)
    sat_pressure_t=1.D0/sat_pressure_t
!    if(SSATG_LOC(ng)<eps)then  ! pure liquid phase
!       if(p_vapor<sat_pressure)then
!          qu_rate=-grid%qu_kin*(p_vapor-sat_pressure)
!          tempvar(1)= -grid%qu_kin !* xxgw
!           tempvar(2)= grid%qu_kin* sat_pressure_t
!          tempvar(3)= 0.D0 ! grid%qu_kin*PPRESSURE_LOC(ng)
!          tempvar(4)=0.D0
!       endif 
!    else if( SSATG_LOC(ng)<(1.D0-eps))then            !2 phases
 !      qu_rate=-grid%qu_kin*(p_vapor-sat_pressure)
!    if( SSATG_LOC(ng)<(1.D0-eps).and. SSATG_LOC(ng)>eps )then
         tempvar(1)= -grid%qu_kin * (xxgw -grid%yh2o_in_co2)
         tempvar(2)= 0.D0!grid%qu_kin* sat_pressure_t
         tempvar(3)= grid%qu_kin*PPRESSURE_LOC(ng)
         tempvar(4)=0.D0

!    else                         !pure gas phase
!          if(p_vapor>sat_pressure .and. SSATG_LOC(ng)>(1.D0-eps))then   
!             tempvar(1)= -grid%qu_kin * xxgw
!             tempvar(2)= grid%qu_kin* sat_pressure_t
!             tempvar(3)=grid%qu_kin*PPRESSURE_LOC(ng)
!             tempvar(4)=0.D0
!          endif
!        end if

 !  if(SSATG_LOC(ng)>eps)then
  !    print *,'n,sg ::',n,SSATG_LOC(ng)
      ra(3,1:4) = ra(3,1:4)-tempvar(1:4)
!   else
!      ra(5,1:4) = ra(5,1:4)+tempvar(1:4)
!      ra(:,4)=0.D0
!      ra(3,1:4)=0.D0
!     ra(3,4)=1.D0
   end if
  
   tempvar=0.D0
   tempvar(1:4)=ra(3,1:4)
   ra(3,1:4)=ra(4,1:4)
   ra(4,1:4)=tempvar(1:4)
   
!  do ii=1,4
!    print *,'accum r',ii,(jj,ra(ii,jj),jj=1,8)
!  enddo
!  end if
   if (grid%iblkfmt == 0) then
     do ii=0,3
      do jj=0,3
        call MatSetValuesLocal(A,1,p1+ii,1,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
      enddo
     enddo
   else
     blkmat11=ra(1:4,1:4)
     call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat11,ADD_VALUES,ierr)
   endif
!print *,'accum r',ra(1:4,1:8)
end do


! -----------------------------contribution from transport----------------------

 !print *,'phase cond: ',iphase_loc_p
  do nc = 1, grid%nconn  ! For each interior connection...
    ra=0.D0
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping   

 
    p1 =  (m1-1)*grid%ndof; t1 = p1+1; c1 = t1+1; s1=c1 + grid%npricomp
    p2 =  (m2-1)*grid%ndof; t2 = p2+1; c2 = t2+1; s2=c2 + grid%npricomp
   
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)  ! determine the normal direction of interface 
    ip2 = grid%iperm2(nc)
    
    iiphas1 = iphase_loc_p(m1)
    iiphas2 = iphase_loc_p(m2)

    select case(ip1)
    case(1) 
       perm1 = perm_xx_loc_p(m1)
    case(2)
       perm1 = perm_yy_loc_p(m1)
    case(3)
       perm1 = perm_zz_loc_p(m1)
    end select
    
    select case(ip2)
    case(1) 
       perm2 = perm_xx_loc_p(m2)
    case(2)
       perm2 = perm_yy_loc_p(m2)
    case(3)
       perm2 = perm_zz_loc_p(m2)
    end select


    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    fluxl = 0.d0

  
! Calculate Henry Coefficient for N2
    fluxlh = 0.d0
    fluxlv = 0.d0
    fluxg = 0.d0
    fluxgh = 0.d0
    fluxgv = 0.d0

    D1 = perm1 !* viscosity_loc_p(m2) 
    D2 = perm2 !* viscosity_loc_p(m1) 
    den = dd2*D1 + dd1*D2
    Dq = (perm1 * perm2) / den   ! average k

 
  ! contribution from rock heat conduction
              i1 = ithrm_loc_p(m1)
              i2 = ithrm_loc_p(m2)
              D1 = grid%ckwet(i1)
              D2 = grid%ckwet(i2)
              Dk = (D1 * D2) / (dd2*D1 + dd1*D2)
              cond = Dk * grid%area(nc)
             
              ra(2,2)= ra(2,2)+cond
              ra(2,6)= ra(2,6)-cond
!print *, 'In 2 ph Jacobian ::  Finished Solid conduction',ra(2,:)

! contribution from diffusion
       por1 = porosity_loc_p(m1)
       por2 = porosity_loc_p(m2)
  
    !in liquid phase
     if(((1.D0-SSATG_LOC(m1))>eps).and.((1.D0-SSATG_LOC(m2))>eps))then

        difff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
        difff = difff * grid%area(nc)*0.25d0
        j=1
        jm1=j+ (m1-1) * grid%nphase
        jm2=j+ (m2-1) * grid%nphase
       
        i1_hencoeff=2+(j-1)*grid%nspec+(m1-1)*grid%nphase*grid%nspec
        i2_hencoeff=2+(j-1)*grid%nspec+(m2-1)*grid%nphase*grid%nspec
        i1_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(m1-1)*grid%nphase*grid%nspec*&
                    grid%npricomp ! +(jc-1)*grid%npricomp
        i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(m2-1)*grid%nphase*grid%nspec*&
                    grid%npricomp ! +(jc-1)*grid%npricomp

       cw1 =1.D0-(CCONC_LOC(m1))*Hen_loc_p(i1_hencoeff)
       cw2 =1.D0-(CCONC_LOC(m2))*Hen_loc_p(i2_hencoeff)
       
       cw1_p= -(CCONC_LOC(m1))*Hen_p_loc_p(i1_hencoeff)
       cw2_p= -(CCONC_LOC(m2))*Hen_p_loc_p(i2_hencoeff)
       cw1_t= -(CCONC_LOC(m1))*Hen_t_loc_p(i1_hencoeff)
       cw2_t= -(CCONC_LOC(m2))*Hen_t_loc_p(i2_hencoeff)
       cw1_c= -Hen_loc_p(i1_hencoeff)-(CCONC_LOC(m1))*Hen_c_loc_p(i1_hencoeff_dc)
       cw2_c= -Hen_loc_p(i2_hencoeff)-(CCONC_LOC(m2))*Hen_c_loc_p(i2_hencoeff_dc)
       cw1_s= -(CCONC_LOC(m1))*Hen_s_loc_p(i1_hencoeff)
       cw2_s= -(CCONC_LOC(m2))*Hen_s_loc_p(i2_hencoeff)


        density_ave = ddensity_loc_p(jm1)+ddensity_loc_p(jm2)
        sat_ave=2.d0 - SSATG_LOC(m1) - SSATG_LOC(m2)
        diffl=difff * sat_ave * density_ave
           

!      difluxl = diffl * grid%area(nc) * (CCONC_LOC(m2) - CCONC_LOC(m1))
        tempvar=0.
        tempvar(1)=difff*sat_ave* &
                ((cw2 - cw1)*d_p_loc_p(jm1)- density_ave *cw1_p)
        tempvar(5)=difff*sat_ave* &
                ((cw2 - cw1)*d_p_loc_p(jm2) +density_ave *cw2_p)
        tempvar(2)=difff*sat_ave* &
                ((cw2 - cw1)*d_t_loc_p(jm1)-density_ave *cw1_t)
        tempvar(6)=difff*sat_ave* &
                ((cw2 - cw1)*d_t_loc_p(jm2) +density_ave *cw2_t)
        tempvar(3)=  difff*sat_ave* &
                ((cw2-cw1)*d_c_loc_p(jm1)-density_ave *cw1_c)
        tempvar(7)=  difff*sat_ave* &
                ((cw2 - cw1)*d_c_loc_p(jm2)+density_ave *cw2_c)
        tempvar(4)=difff *(- (cw2 - cw1) *density_ave + d_s_loc_p(jm1)* &
                sat_ave*(cw2 - cw1)-cw1_s * sat_ave*density_ave)
        tempvar(8)=difff *(- (cw2 - cw1) *density_ave + d_s_loc_p(jm2)* &
                sat_ave*(cw2 - cw1)+cw2_s * sat_ave*density_ave)
        
        
!        ra(5,:)= ra(5,:) - tempvar
        ra(4,:)= ra(4,:) - tempvar
 !    print *, 'After diff ', tempvar   
     endif

  !print *,'In 2 ph Jacobian ::  Finished Liquid diff'
   !in gas phase
    if((SSATG_LOC(m1)>eps).and.(SSATG_LOC(m2))>eps)then 
       j=2
        jm1=j+ (m1-1) * grid%nphase
        jm2=j+ (m2-1) * grid%nphase
        difff = (por1 * por2) / (dd2*por1 + dd1*por2) 
        difff = difff * grid%area(nc)*0.25
        diffg=(grid%cdiff(int(ithrm_loc_p(m1)))+grid%cdiff(int(ithrm_loc_p(m2))))*0.5D0
       
        difff = difff *diffg
        density_ave = ddensity_loc_p(jm1)+ddensity_loc_p(jm2)
        sat_ave= SSATG_LOC(m1) + SSATG_LOC(m2)
        
        diffg=difff * density_ave * sat_ave
      !  print *,'Joca, difff, diffg=',difff,diffg, sat_ave,density_ave
      
        uconc=-CCONC_LOC(m2)+ CCONC_LOC(m1)
        tempvar=0.d0
        tempvar(1)=difff*sat_ave * uconc * d_p_loc_p(jm1)
        tempvar(5)=difff*sat_ave* uconc *d_p_loc_p(jm2)
        tempvar(2)=difff*sat_ave* uconc *d_t_loc_p(jm1)
        tempvar(6)=difff*sat_ave* uconc *d_t_loc_p(jm2)
        tempvar(3)=difff*sat_ave*(uconc * &
                   d_c_loc_p(jm1)+density_ave)
        tempvar(7)= difff*sat_ave * (uconc * &
                   d_c_loc_p(jm2)-density_ave)
        tempvar(4)=difff *(density_ave* uconc +&
                   sat_ave * d_s_loc_p(jm1)*uconc )
        tempvar(8)=difff *(density_ave* uconc + &
                   sat_ave * d_s_loc_p(jm2)*uconc )

        ra(3,:)= ra(3,:) - tempvar(:)
        ra(4,:)= ra(4,:) - tempvar(:)

     endif
     
!print *,'In 2 ph Jacobian ::  Finished gas diff'
!print *,'diff ',ra(1:5,1:8) 

!********for liquid flux************************** 
    if(((1.D0-SSATG_LOC(m1))>grid%swir(int(icap_loc_p(m1)))).or.&
         ((1.D0-SSATG_LOC(m2))>grid%swir(int(icap_loc_p(m2)))))then
       j=1             ! index to handle varibles defined with context da_nphase
    ! j=1 means liquid phase
       jm1= j + (m1-1) * grid%nphase
       jm2= j + (m2-1) * grid%nphase


       upweight=0.5d0  ! weight for m1
       if((1.D0-SSATG_LOC(m1)) < grid%swir(int(icap_loc_p(m1))))then
          upweight=0.d0
       else if((1.D0-SSATG_LOC(m2)) < grid%swir(int(icap_loc_p(m2))))then
          upweight=1.d0
       endif

       density_ave = upweight*ddensity_loc_p(jm1)+ &
                     (1.D0-upweight)*ddensity_loc_p(jm2)  

       gravity = (upweight*ddensity_loc_p(jm1)*avgmw_loc_p(jm1) + &
            (1.D0-upweight)*ddensity_loc_p(jm2)*avgmw_loc_p(jm2)) &
            * grid%gravity * grid%delz(nc)
 

       dphi= -(PPRESSURE_LOC(m2) - PPRESSURE_LOC(m1) & 
                 - pc_loc_p(jm2) + pc_loc_p(jm1)   &
                 - gravity)

       if(dphi>=0.D0)then
         mu=m1
       else
         mu=m2
       end if
       jmu = j + (mu-1) * grid%nphase          
       

       if((kvr_loc_p(jmu)*Dq) > floweps)then 

          D0=kvr_loc_p(jmu) 

          v_darcy = Dq * D0 * dphi
          q = v_darcy * grid%area(nc)
          !q = 0.d0
          gdz=grid%gravity * grid%delz(nc)


         
          if(m1==mu)then   
             m1weight = 1.D0
             m2weight = 0.D0
          else
             m1weight = 0.D0
             m2weight = 1.D0
          endif

          devq=0.D0    
          devq(1,1) = Dq*(kvr_p_loc_p(jm1)*dphi*m1weight + &
                       D0* (1.D0-pc_p_loc_p(jm1))+gdz* &
                       upweight*d_p_loc_p(jm1)*avgmw_loc_p(jm1))
          devq(1,2) = Dq*(kvr_p_loc_p(jm2)*dphi*m2weight + &
                       D0* (-1.+pc_p_loc_p(jm2))+gdz* (1.D0- &
                       upweight)*d_p_loc_p(jm2)*avgmw_loc_p(jm2))
           ! volume rate derivative to p1,p2
           
          devq(2,1)=Dq*(kvr_t_loc_p(jm1)*dphi*m1weight + D0* &
                     (-pc_t_loc_p(jm1)+ gdz*upweight*avgmw_loc_p(jm1)&
                     *  d_t_loc_p(jm1)))
          devq(2,2)=Dq*(kvr_t_loc_p(jm2)*dphi*m2weight + D0* &
                     (pc_t_loc_p(jm2)+ gdz*(1.D0-upweight)*   &
                     avgmw_loc_p(jm2) *  d_t_loc_p(jm2)))
           ! volume rate derivative to T1,T2

          devq(3,1)=Dq*(kvr_c_loc_p(jm1)*dphi*m1weight &
                     +D0*(-pc_c_loc_p(jm1)+ gdz*upweight* &
                     (avgmw_loc_p(jm1) *d_c_loc_p(jm1)+ &
                     ddensity_loc_p(jm1)*avgmw_c_loc_p(jm1))))
          devq(3,2)=Dq*(kvr_c_loc_p(jm2)*dphi*m2weight &
                     +D0*(pc_c_loc_p(jm2)+ gdz*(1.D0-upweight)* &
                     (avgmw_loc_p(jm2) *d_c_loc_p(jm2)+ &
                     ddensity_loc_p(jm2)*avgmw_c_loc_p(jm2))))
           ! volume rate derivative to C1,C2

          devq(4,1)=Dq*(kvr_s_loc_p(jm1)*dphi*m1weight - &
                     D0*pc_s_loc_p(jm1))
          devq(4,2)=Dq*(kvr_s_loc_p(jm2)*dphi*m2weight + &
                     D0*pc_s_loc_p(jm2))
           ! volume rate derivative to s1,s2
          devq=devq*grid%area(nc)




! pressure equation
          ra(1,1)=ra(1,1)+devq(1,1)*density_ave+q*upweight*d_p_loc_p(jm1)
          ra(1,5)=ra(1,5)+devq(1,2)*density_ave+q*(1.D0-upweight)*d_p_loc_p(jm2)
          ra(1,2)=ra(1,2)+devq(2,1)*density_ave+q*upweight*d_t_loc_p(jm1)
          ra(1,6)=ra(1,6)+devq(2,2)*density_ave+q*(1.D0-upweight)*d_t_loc_p(jm2)
          ra(1,3)=ra(1,3)+devq(3,1)*density_ave+q*upweight*d_c_loc_p(jm1)
          ra(1,7)=ra(1,7)+devq(3,2)*density_ave+q*(1.D0-upweight)*d_c_loc_p(jm2)
          ra(1,4)=ra(1,4)+devq(4,1)*density_ave+q*upweight*d_s_loc_p(jm1)
          ra(1,8)=ra(1,8)+devq(4,2)*density_ave+q*(1.D0-upweight)*d_s_loc_p(jm2)
       !   print *,' Jacobian:: ra(1,:)  ',density_ave, d_p_loc_p(jm1:jm2)
       !   print *,ra(1,:)

! energy equation
!  remember 
           ra(2,1)= ra(2,1)+hh_loc_p(jmu)*density_ave*devq(1,1) &
                   + q * hh_loc_p(jmu)*d_p_loc_p(jm1)*upweight &
                   + q * density_ave * h_p_loc_p(jm1)* m1weight
           ra(2,5)= ra(2,5)+hh_loc_p(jmu)*density_ave*devq(1,2) &
                   + q * hh_loc_p(jmu)*d_p_loc_p(jm2)*(1.D0-upweight) &
                   + q * density_ave * h_p_loc_p(jm2)* m2weight
           ra(2,2)= ra(2,2)+hh_loc_p(jmu)*density_ave*devq(2,1) &
                   + q * hh_loc_p(jmu)*d_t_loc_p(jm1)*upweight &
                   + q * density_ave * h_t_loc_p(jm1)* m1weight
           ra(2,6)= ra(2,6)+ hh_loc_p(jmu)*density_ave*devq(2,2) &
                   + q * hh_loc_p(jmu)*d_t_loc_p(jm2)*(1.D0-upweight) &
                   + q * density_ave * h_t_loc_p(jm2)* m2weight
           ra(2,3)= ra(2,3)+hh_loc_p(jmu)*density_ave*devq(3,1) &
                   + q * hh_loc_p(jmu)*d_c_loc_p(jm1)*upweight &
                   + q * density_ave * h_c_loc_p(jm1)* m1weight
           ra(2,7)= ra(2,7)+hh_loc_p(jmu)*density_ave*devq(3,2) &
                   + q * hh_loc_p(jmu)*d_c_loc_p(jm2)*(1.D0-upweight) &
                   + q * density_ave * h_c_loc_p(jm2)* m2weight
           ra(2,4)= ra(2,4)+hh_loc_p(jmu)*density_ave*devq(4,1) &
                   + q * hh_loc_p(jmu)*d_s_loc_p(jm1)*upweight &
                   + q * density_ave * h_s_loc_p(jm1)* m1weight
           ra(2,8)= ra(2,8)+hh_loc_p(jmu)*density_ave*devq(4,2) &
                   + q * hh_loc_p(jmu)*d_s_loc_p(jm2)*(1.D0-upweight) &
                   + q * density_ave * h_s_loc_p(jm2)* m2weight
 
 
! Saturation Equation 
           iu_hencoeff=2+(j-1)*grid%nspec+(mu-1)*grid%nphase*grid%nspec
           iu_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(mu-1)*grid%nphase*&
                          grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp

           i1_hencoeff=2+(j-1)*grid%nspec+(m1-1)*grid%nphase*grid%nspec
           i1_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(m1-1)*grid%nphase*&
                          grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp

 
! Calculate Henry Coefficient for N2
           i2_hencoeff=2+(j-1)*grid%nspec+(m2-1)*grid%nphase*grid%nspec
           i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(m2-1)*grid%nphase*&
                          grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp

!   fluxgv = fluxgv + density_ave * q * (Henny_loc_p(1+(j-1)*grid%nspec+ &
!                (jmu-1)*grid%nphase*grid%nspec)*(1-CCONC_LOC(mu)))
           
          
           cw1 =1.D0-(CCONC_LOC(m1))*Hen_loc_p(i1_hencoeff)
           cw2 =1.D0-(CCONC_LOC(m2))*Hen_loc_p(i2_hencoeff)
       
           cw1_p= -(CCONC_LOC(m1))*Hen_p_loc_p(i1_hencoeff)
           cw2_p= -(CCONC_LOC(m2))*Hen_p_loc_p(i2_hencoeff)
           cw1_t= -(CCONC_LOC(m1))*Hen_t_loc_p(i1_hencoeff)
           cw2_t= -(CCONC_LOC(m2))*Hen_t_loc_p(i2_hencoeff)
           cw1_c= -Hen_loc_p(i1_hencoeff)-(CCONC_LOC(m1))*Hen_c_loc_p(i1_hencoeff_dc)
           cw2_c= -Hen_loc_p(i2_hencoeff)-(CCONC_LOC(m2))*Hen_c_loc_p(i2_hencoeff_dc)
           cw1_s= -(CCONC_LOC(m1))*Hen_s_loc_p(i1_hencoeff)
           cw2_s= -(CCONC_LOC(m2))*Hen_s_loc_p(i2_hencoeff)

           if(m1==mu)then
              cwu=cw1
             else
               cwu=cw2  
            end if

           tempvar = 0.D0 
           tempvar(1) = d_p_loc_p(jm1)*upweight *q *cwu &
                + density_ave *devq(1,1) * cwu   &
                + density_ave * q * m1weight * cw1_p
           tempvar(5) = d_p_loc_p(jm2)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(1,2) * cwu   &
                + density_ave * q * m2weight *cw2_p

           tempvar(2) = d_t_loc_p(jm1)*upweight *q *cwu &
                + density_ave *devq(2,1) * cwu   &
                + density_ave * q * m1weight * cw1_t
           tempvar(6) = d_t_loc_p(jm2)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(2,2) * cwu   &
                + density_ave * q * m2weight * cw2_t

           tempvar(3) = d_c_loc_p(jm1)*upweight *q * cwu &
                + density_ave *devq(3,1) * cwu   &
                + density_ave * q * m1weight * cw1_c
           tempvar(7)= d_c_loc_p(jm2)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(3,2) * cwu   &
                + density_ave * q * m2weight*cw2_c

           tempvar(4) = d_s_loc_p(jm1)*upweight *q * cwu &
                + density_ave *devq(4,1) * cwu   &
                + density_ave * q * m1weight * cw1_s
           tempvar(8) = d_s_loc_p(jm2)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(4,2) * cwu   &
                + density_ave * q * m2weight * cw2_s

           ra(4,:) = ra(4,:) + tempvar(:)
 !          ra(5,:) = ra(5,:) + tempvar(:)

        end if
     end if

!print *,'In 2ph Jacobian :: Finished Liquid Convection ',ra(4,:)
!************* for Gas phase Flux************************************************* 


     if((SSATG_LOC(m1)>eps).or.(SSATG_LOC(m2)>eps))then

        j=2             ! index to handle varibles defined with context da_nphase
        jm1= j + (m1-1) * grid%nphase
        jm2= j + (m2-1) * grid%nphase
                  
        upweight=0.5D0  ! weight for m1
        if(SSATG_LOC(m1)<eps)then
           upweight=0.D0
        else if(SSATG_LOC(m2)<eps)then
           upweight=1.D0
        endif
        density_ave = upweight*ddensity_loc_p(jm1)+ &
             (1.D0-upweight)*ddensity_loc_p(jm2)  

        gravity = (upweight*ddensity_loc_p(jm1)*avgmw_loc_p(jm1) + &
              (1.D0-upweight)*ddensity_loc_p(jm2)*avgmw_loc_p(jm2)) &
              * grid%gravity * grid%delz(nc)

        dphi= -PPRESSURE_LOC(m2) + PPRESSURE_LOC(m1)+ gravity
! pressure gradient does not contain pc term

        if(dphi>=0.D0)then
           mu=m1
        else
           mu=m2
        end if
        jmu= j + (mu-1) * grid%nphase        
       ! print *,'Jac, Gas q::',m1,m2,dphi,mu
        !print *, ' Gas flow mobility::', kvr_loc_p(jmu)*Dq
        if((kvr_loc_p(jmu)*Dq)>floweps)then
           D0=kvr_loc_p(jmu)
           v_darcy =  Dq * D0 * dphi
           gdz=grid%gravity * grid%delz(nc)
           q = v_darcy * grid%area(nc)
           !q = 0.d0

           if(m1==mu)then   
              m1weight = 1.D0
              m2weight = 0.D0
           else
              m1weight = 0.D0
              m2weight = 1.D0
           endif

           devq(1,1) = Dq*(kvr_p_loc_p(jm1)*dphi*m1weight + &
                D0* (1.D0-pc_p_loc_p(jm1))+gdz* &
                upweight*d_p_loc_p(jm1)*avgmw_loc_p(jm1))
           devq(1,2) = Dq*(kvr_p_loc_p(jm2)*dphi*m2weight + &
                D0* (-1.d0+pc_p_loc_p(jm2))+gdz* (1.D0- &
                upweight)*d_p_loc_p(jm2)*avgmw_loc_p(jm2))
           ! volume rate derivative to p1,p2
           
           devq(2,1)=Dq*(kvr_t_loc_p(jm1)*dphi*m1weight + D0* &
                (-pc_t_loc_p(jm1)+ gdz*upweight*avgmw_loc_p(jm1)&
                *  d_t_loc_p(jm1)))
           devq(2,2)=Dq*(kvr_t_loc_p(jm2)*dphi*m2weight + D0* &
                (pc_t_loc_p(jm2)+ gdz*(1.D0-upweight)*   &
                avgmw_loc_p(jm2) *  d_t_loc_p(jm2)))
           ! volume rate derivative to T1,T2

           devq(3,1)=Dq*(kvr_c_loc_p(jm1)*dphi*m1weight &
                +D0*(-pc_c_loc_p(jm1)+ gdz*upweight* &
                (avgmw_loc_p(jm1) *d_c_loc_p(jm1)+ &
                ddensity_loc_p(jm1)*avgmw_c_loc_p(jm1))))
           devq(3,2)=Dq*(kvr_c_loc_p(jm2)*dphi*m2weight &
                +D0*(pc_c_loc_p(jm2)+ gdz*(1.D0-upweight)* &
                (avgmw_loc_p(jm2) *d_c_loc_p(jm2)+ &
                ddensity_loc_p(jm2)*avgmw_c_loc_p(jm2))))
           ! volume rate derivative to C1,C2

           devq(4,1)=Dq*(kvr_s_loc_p(jm1)*dphi*m1weight - &
                D0*pc_s_loc_p(jm1))
           devq(4,2)=Dq*(kvr_s_loc_p(jm2)*dphi*m2weight + &
                D0*pc_s_loc_p(jm2))

           devq=devq* grid%area(nc)
 !print *,'devq:',nc,q,dphi,devq(3,:)
         !  print *,'Jac devg',devq(3,:)
! volume rate derivative to s1,s2
! the only difference between gas and water phase is the pc terms
! by setting pc=0 and dpc=0 for gas phase, they are all the same

 
! Calculate Henry Coefficient for N2
           !print *,' Jacobian::dq',devq
    !       print *,'ra  orig', ra(1,:)
           ra(1,1)=ra(1,1)+devq(1,1)*density_ave+q*upweight*d_p_loc_p(jm1)
           ra(1,5)=ra(1,5)+devq(1,2)*density_ave+q*(1.D0-upweight)*d_p_loc_p(jm2)
           ra(1,2)=ra(1,2)+devq(2,1)*density_ave+q*upweight*d_t_loc_p(jm1)
           ra(1,6)=ra(1,6)+devq(2,2)*density_ave+q*(1.D0-upweight)*d_t_loc_p(jm2)
           ra(1,3)=ra(1,3)+devq(3,1)*density_ave+q*upweight*d_c_loc_p(jm1)
           ra(1,7)=ra(1,7)+devq(3,2)*density_ave+q*(1.D0-upweight)*d_c_loc_p(jm2)
           ra(1,4)=ra(1,4)+devq(4,1)*density_ave+q*upweight*d_s_loc_p(jm1)
           ra(1,8)=ra(1,8)+devq(4,2)*density_ave+q*(1.D0-upweight)*d_s_loc_p(jm2)
     !      print *, 'jm ', jm1, jm2
      !     print *,' Jacobian:: den  ',density_ave, d_p_loc_p(jm1),d_p_loc_p(jm2)
      !     print *,' Jacobian:: den  ', d_t_loc_p(jm1), d_t_loc_p(jm2), d_c_loc_p(jm1),d_c_loc_p(jm2)
           !print *,' Jacobian:: ra(1,:)  ',kvr_c_loc_p(jm1),kvr_c_loc_p(jm2),pc_c_loc_p(jm1),pc_c_loc_p(jm2)
           !print *,' Jacobian:: ra(1,:)q=  ',q,gdz
         !  print *,'devq:',nc,q,dphi,devq(1,:)
!          print *,ra(1,:)
! terms pressure equation the same for every phase
    ! print *,'RA(2,3;7) orig: ',ra(2,3) ,ra(2,7)
    ! print *, 'Devq ',devq(3,1),devq(3,2),q
    ! print *, 'Weight ',upweight, m1weight, m2weight,jm1,jm2,jmu
    ! print *, 'Dev_c ',d_c_loc_p(jm1),d_c_loc_p(jm2),h_c_loc_p(jm1),h_c_loc_p(jm2)
           ra(2,1)=ra(2,1)+ hh_loc_p(jmu)*density_ave*devq(1,1) &
                + q * hh_loc_p(jmu)*d_p_loc_p(jm1)*upweight &
                + q * density_ave * h_p_loc_p(jm1)* m1weight
           ra(2,5)=ra(2,5)+ hh_loc_p(jmu)*density_ave*devq(1,2) &
                + q * hh_loc_p(jmu)*d_p_loc_p(jm2)*(1.D0-upweight) &
                + q * density_ave * h_p_loc_p(jm2)* m2weight
           ra(2,2)=ra(2,2)+ hh_loc_p(jmu)*density_ave*devq(2,1) &
                + q * hh_loc_p(jmu)*d_t_loc_p(jm1)*upweight &
                + q * density_ave * h_t_loc_p(jm1)* m1weight
           ra(2,6)=ra(2,6)+ hh_loc_p(jmu)*density_ave*devq(2,2) &
                + q * hh_loc_p(jmu)*d_t_loc_p(jm2)*(1.D0-upweight) &
                + q * density_ave * h_t_loc_p(jm2)* m2weight
           ra(2,3)=ra(2,3)+ hh_loc_p(jmu)*density_ave*devq(3,1) &
                + q * hh_loc_p(jmu)*d_c_loc_p(jm1)*upweight &
                + q * density_ave * h_c_loc_p(jm1)* m1weight
           ra(2,7)=ra(2,7)+ hh_loc_p(jmu)*density_ave*devq(3,2) &
                + q * hh_loc_p(jmu)*d_c_loc_p(jm2)*(1.D0-upweight) &
                + q * density_ave * h_c_loc_p(jm2)* m2weight
           ra(2,4)=ra(2,4)+ hh_loc_p(jmu)*density_ave*devq(4,1) &
                + q * hh_loc_p(jmu)*d_s_loc_p(jm1)*upweight &
                + q * density_ave * h_s_loc_p(jm1)* m1weight
           ra(2,8)=ra(2,8)+ hh_loc_p(jmu)*density_ave*devq(4,2) &
                + q * hh_loc_p(jmu)*d_s_loc_p(jm2)*(1.D0-upweight) &
                + q * density_ave * h_s_loc_p(jm2)* m2weight

   !print *,'RA(2,3) now 1: ',ra(2,3),ra(2,7)

! terms in energy equation the same for every phase

! concentration equation
           tempvar=0.
           uconc = 1.D0- CCONC_LOC(mu)
           tempvar(1)= uconc*density_ave*devq(1,1) &
                   + q * uconc*d_p_loc_p(jm1)*upweight
                   
           tempvar(5)= uconc*density_ave*devq(1,2) &
                   + q * uconc*d_p_loc_p(jm2)*(1.D0-upweight)
                  
           tempvar(2)= uconc*density_ave*devq(2,1) &
                   + q * uconc*d_t_loc_p(jm1)*upweight
           tempvar(6)= uconc*density_ave*devq(2,2) &
                   + q * uconc*d_t_loc_p(jm2)*(1.D0-upweight)

           tempvar(3)= uconc*density_ave*devq(3,1)&
                   + q * uconc*d_c_loc_p(jm1)*upweight &
                   - q * density_ave * m1weight
           tempvar(7)= uconc*density_ave*devq(3,2) &
                   + q * uconc*d_c_loc_p(jm2)*(1.D0-upweight) &
                   - q * density_ave * m2weight    
           
           tempvar(4)= uconc*density_ave*devq(4,1) &
                   + q * uconc*d_s_loc_p(jm1)*upweight
           tempvar(8)= uconc*density_ave*devq(4,2) &
                   + q * uconc*d_s_loc_p(jm2)*(1.D0-upweight)

           ra(3,:)=ra(3,:)+tempvar(:)
          
! Saturation Equation 
           ra(4,:)=ra(4,:)+tempvar(:)
!           print *,'Gas ',uconc,q,devq(1,1), devq(1,2),tempvar
        end if
     end if
     !print *,'In 2ph Jacobian :: Finished Gas Convection '

!     if(SSATG_LOC(m1)<eps) ra(:,4)=0.D0
!     if(SSATG_LOC(m2)<eps) ra(:,8)=0.D0

tempvar=0.D0

     tempvar(1:8)=ra(3,1:8)
     ra(3,1:8)=ra(4,1:8)
     ra(4,1:8)=tempvar(1:8)

 ! if(SSATG_LOC(ng)<=eps)then
 !    ra(4,1:8)=ra(5,1:8)
 ! endif   
   
  ! print *,'Flux r 90: ',ra(1:4,1:8)
  ! do ii=1,4
  !   print *,'Flux r',ii,(jj,':',ra(ii,jj),jj=1,8)
  ! enddo
  
   if(grid%iblkfmt == 1) then
     blkmat11 = 0.D0; blkmat12 = 0.D0; blkmat21 = 0.D0; blkmat22 = 0.D0;
   endif
   
     do ii=0,3
       do jj=0,3
         select case(ii)
         case (3)
           if(n1>0)then 
              if(iiphas1==6)then
                if (grid%iblkfmt == 0) then
                  call MatSetValuesLocal(A,1,p1+ii,1,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
                else
                  blkmat11(ii+1,jj+1) =blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
                endif
     !        else
     !          call  MatSetValuesLocal(A,1,p1+ii,1,p1+jj,ra(5,jj+1),ADD_VALUES,ierr)
              endif
           end if
           if(n2>0)then 
             if(iiphas2==6)then
                if (grid%iblkfmt == 0) then
                  call MatSetValuesLocal(A,1,p2+ii,1,p1+jj,-ra(ii+1,jj+1),ADD_VALUES,ierr)
                else
                  blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1)-ra(ii+1,jj+1)
                endif
          !  else
         !     call MatSetValuesLocal(A,1,p2+ii,1,p1+jj,-ra(5,jj+1),ADD_VALUES,ierr)
             end if
          end if
        case default      
           if(n1>0) then
             if (grid%iblkfmt == 0) then
               call MatSetValuesLocal(A,1,p1+ii,1,p1+jj,ra(ii+1,jj+1),ADD_VALUES,ierr)
             else
               blkmat11(ii+1,jj+1) = blkmat11(ii+1,jj+1) + ra(ii+1,jj+1)
             endif
           endif
           if(n2>0) then
             if (grid%iblkfmt == 0) then
               call MatSetValuesLocal(A,1,p2+ii,1,p1+jj,-ra(ii+1,jj+1),ADD_VALUES,ierr)
             else
               blkmat21(ii+1,jj+1) = blkmat21(ii+1,jj+1) -ra(ii+1,jj+1)
             endif
           endif
        end select
        enddo
        do jj=4,7
           select case(ii)
           case (3)
              if(n1>0)then 
                if(iiphas1==6)then
                  if (grid%iblkfmt == 0) then
                    call MatSetValuesLocal(A,1,p1+ii,1,p2+jj-4,ra(ii+1,jj+1),ADD_VALUES,ierr)
                  else
                    blkmat12(ii+1,jj-4+1) = blkmat12(ii+1,jj-4+1) + ra(ii+1,jj+1)
                  endif
            !    else
            !      call MatSetValuesLocal(A,1,p1+ii,1,p2+jj-4,ra(5,jj+1),ADD_VALUES,ierr)
                 end if
              end if
              if(n2>0)then 
                if(iiphas2==6)then
                  if (grid%iblkfmt == 0) then
                    call MatSetValuesLocal(A,1,p2+ii,1,p2+jj-4,-ra(ii+1,jj+1),ADD_VALUES,ierr)
                  else
                    blkmat22(ii+1,jj-4+1) =  blkmat22(ii+1,jj-4+1) - ra(ii+1,jj+1)
                  endif
             !   else
             !     call MatSetValuesLocal(A,1,p2+ii,1,p2+jj-4,-ra(5,jj+1),ADD_VALUES,ierr)
                 end if
              end if
           case default   
              if(n1>0) then
                  if (grid%iblkfmt == 0) then
                    call MatSetValuesLocal(A,1,p1+ii,1,p2+jj-4,ra(ii+1,jj+1),ADD_VALUES,ierr)
                  else
                    blkmat12(ii+1,jj-4+1) = blkmat12(ii+1,jj-4+1) + ra(ii+1,jj+1)
                  endif
              endif
              if(n2>0) then
                  if (grid%iblkfmt == 0) then
                    call MatSetValuesLocal(A,1,p2+ii,1,p2+jj-4,-ra(ii+1,jj+1),ADD_VALUES,ierr)
                  else
                    blkmat22(ii+1,jj-4+1) =  blkmat22(ii+1,jj-4+1) - ra(ii+1,jj+1)
                  endif
              endif
           end select
        enddo
    enddo
  
  if (grid%iblkfmt /= 0) then
    call MatSetValuesBlockedLocal(A,1,m1-1,1,m1-1,blkmat11,ADD_VALUES,ierr)
    call MatSetValuesBlockedLocal(A,1,m2-1,1,m2-1,blkmat22,ADD_VALUES,ierr)
    call MatSetValuesBlockedLocal(A,1,m1-1,1,m2-1,blkmat12,ADD_VALUES,ierr)
    call MatSetValuesBlockedLocal(A,1,m2-1,1,m1-1,blkmat21,ADD_VALUES,ierr)
  endif
!print *,'accum r',ra(1:5,1:8)   
 !print *,'devq:',nc,q,dphi,devq(3,:)
  end do

!************** handle boundary conditions ***************************8

  do nc = 1, grid%nconnbc

     ra=0.
     m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
     ng = grid%nL2G(m)
     
     p1 = (ng-1) * grid%ndof
     t1 = p1 + 1
     c1 = t1 + 1
     s1=  c1 + 1
     ibc = grid%ibconn(nc)
     ip1 = grid%ipermbc(nc) 
     iiphas = iphase_loc_p(ng)

     select case(ip1)
     case(1)
       perm1 = perm_xx_loc_p(ng)
     case(2)
       perm1 = perm_yy_loc_p(ng)
     case(3)
       perm1 = perm_zz_loc_p(ng)
     end select

    select case(grid%ibndtyp(ibc))
      case(1)
      case(2)
      ! solve for pb from Darcy's law given qb /= 0
        grid%pressurebc(:,nc)=PPRESSURE_LOC(ng)
        grid%tempbc(nc) = TTEMP_LOC(ng)
        grid%sgbc(nc) = SSATG_LOC(ng)
        grid%concbc(nc) = CCONC_LOC(ng)

      case(3) 
        grid%tempbc(nc) = TTEMP_LOC(ng)
        grid%sgbc(nc) = SSATG_LOC(ng)
        grid%concbc(nc) = CCONC_LOC(ng)
    end select

!   print *,'pflow_2pha-bc2: ',ibc,grid%ideriv,grid%density_bc
    
   if (grid%ideriv == 1) then
      iicap = icap_loc_p(ng)     
      call mixture_eos(grid%pressurebc(2,nc),grid%tempbc(nc),&
      grid%concbc(nc),grid%sgbc(nc),&
      grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
      iicap, grid%swir(iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&       !use the node's value
      grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,grid%density_bc,&
      grid%d_p_bc,grid%d_t_bc,grid%d_c_bc,grid%d_s_bc,&
      grid%avgmw_bc,grid%avgmw_c_bc,&
      grid%hh_bc,grid%h_p_bc,grid%h_t_bc,grid%h_c_bc,grid%h_s_bc,&
      grid%uu_bc,grid%u_p_bc,grid%u_t_bc,grid%u_c_bc,grid%u_s_bc,&
      grid%df_bc,grid%df_p_bc,grid%df_t_bc,grid%df_s_bc,grid%df_c_bc,&
      grid%hen_bc,grid%hen_p_bc,grid%hen_t_bc,grid%hen_s_bc,grid%hen_c_bc,&
      grid%xxphi_co2_bc(nc),&
      grid%pc_bc,grid%pc_p_bc,grid%pc_t_bc,grid%pc_c_bc,grid%pc_s_bc,&
      grid%kvr_bc,grid%kvr_p_bc,grid%kvr_t_bc,grid%kvr_c_bc,grid%kvr_s_bc,&
      ierr,grid%itable)
   end if
  
      select case (grid%ibndtyp(ibc))
      case(1)   ! Dirichlet BC for p, T, C, S

 !Cond and diff     
         i1 = ithrm_loc_p(ng)
         cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
         ra(2,6)=ra(2,6)-cond
    !     print *,'BC cond',cond
         por2 = porosity_loc_p(ng)
         
         ! diff in gas  phase
          if((grid%sgbc(nc)>eps).and.(SSATG_LOC(ng)>eps))then 
             j=2
             jng=j+ (ng-1) * grid%nphase
             
             diffg=grid%cdiff(int(ithrm_loc_p(ng)))
             difff = por2  * diffg * grid%areabc(nc) * 0.25/ grid%distbc(nc)
             density_ave=grid%density_bc(j)+ddensity_loc_p(jng)
             sat_ave=SSATG_LOC(ng)+grid%sgbc(nc)
             diffg=difff *sat_ave * density_ave
 
             
             tempvar=0.D0
             acc=(-CCONC_LOC(ng) + grid%concbc(nc))
             tempvar(5)=difff * sat_ave * acc * d_p_loc_p(jng)
             tempvar(6)=difff* sat_ave * acc*d_t_loc_p(jng)
             tempvar(7)= - diffg + difff * sat_ave * acc*d_c_loc_p(jng)
             tempvar(8)= difff *acc * density_ave

             ra(3,5:8)= ra(3,5:8) - tempvar(5:8)
             ra(4,5:8)= ra(4,5:8) - tempvar(5:8)
          end if

          !diff in Liquid  phase
          tempvar=0.D0
          if(((1.D0-grid%sgbc(nc))>eps).and.((1.D0-SSATG_LOC(ng))>eps))then 
             j=1
             jng =j+ (ng-1) * grid%nphase

              diffl=grid%difaq / grid%distbc(nc)
              difff = por2  * diffl * grid%areabc(nc) * 0.25
              density_ave=grid%density_bc(j)+ddensity_loc_p(jng)
              sat_ave=2.D0-SSATG_LOC(ng)-grid%sgbc(nc)
              diffl=difff *sat_ave * density_ave
 
             ibc_hencoeff=2+(j-1)*grid%nspec
!             ibc_hencoeff_dc=1+(j-1)*grid%npricomp*grid%nspec
             i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
             i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(ng-1)*grid%nphase*&
                  grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp



             cw1 =1.D0-(grid%concbc(nc))*grid%Hen_bc(ibc_hencoeff)
             cw2 =1.D0-(CCONC_LOC(ng))*Hen_loc_p(i2_hencoeff)

             cw2_p= -(CCONC_LOC(ng))*Hen_p_loc_p(i2_hencoeff)
             cw2_t= -(CCONC_LOC(ng))*Hen_t_loc_p(i2_hencoeff)
             cw2_c= -Hen_loc_p(i2_hencoeff)-(CCONC_LOC(ng))*Hen_c_loc_p(i2_hencoeff_dc)
             cw2_s= -(CCONC_LOC(ng))*Hen_s_loc_p(i2_hencoeff)

             tempvar(5)=difff*sat_ave* &
                  ((cw2 - cw1)*d_p_loc_p(jng) +density_ave *cw2_p)
             tempvar(6)=difff*sat_ave* &
                  ((cw2 - cw1)*d_t_loc_p(jng) +density_ave *cw2_t)
             tempvar(7)=  difff*sat_ave* &
                  ((cw2 - cw1)*d_c_loc_p(jng)+density_ave *cw2_c)
             tempvar(8)=difff *(- (cw2 - cw1) *density_ave + d_s_loc_p(jng)* &
                  sat_ave*(cw2 - cw1)+cw2_s * sat_ave*density_ave)


             ra(4,5:8)= ra(4,5:8) - tempvar(5:8)
     !        ra(5,5:8)= ra(5,5:8) - tempvar(5:8)
     !        print *,'     In 2ph Jacobian ::  BC 1 : diff and cond term '
     !        print *,grid%concbc(ibc), cw1,cw2, tempvar(5:8)
     !        print *,d_p_loc_p(jng),cw2_p
         end if
 
!print *,ra   

! Convection terms
         Dq= perm1 / grid%distbc(nc)

      ! Liquid phase   
         if(((1.D0-grid%sgbc(nc))>grid%swir(int(icap_loc_p(ng)))).or.    &
                ((1.D0-SSATG_LOC(ng))>grid%swir(int(icap_loc_p(ng)))))then
          
            j=1
            jng = j + (ng-1) * grid%nphase 

            upweight=1.D0
            if((1.D0-grid%sgbc(nc))<grid%swir(int(icap_loc_p(ng))))then
               upweight=0.d0
            endif
            density_ave = upweight* grid%density_bc(j)+ &
                        (1.D0-upweight)*ddensity_loc_p(jng)  

            gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                 (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                 * grid%gravity * grid%delzbc(nc)
            dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(1,nc) & 
                 - pc_loc_p(jng) + grid%pc_bc(j)   &
                 - gravity)
            
          if(dphi>=0.D0)then
             mu=0 ! BC
             ukvr=grid%kvr_bc(j)
             uhh=grid%hh_bc(j)
             uconc=grid%concbc(nc)
            else 
             mu=ng
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc=CCONC_LOC(ng)
          endif

           if((ukvr*Dq)>floweps)then 

              D0=ukvr 

              v_darcy =  Dq * D0 * dphi
              q = v_darcy * grid%areabc(nc)
              !q = 0.d0
              gdz=grid%gravity * grid%delzbc(nc)
                              
              if(mu==0)then   
                 mbweight = 1.D0
                 mnweight = 0.D0
              else
                 mbweight = 0.D0
                 mnweight = 1.D0
              endif

! now 1 refer to boundary
!     2          node contains BC
! So do not take derivative to index 2
          devq=0.D0    
          devq(1,2) = Dq*(kvr_p_loc_p(jng)*dphi*mnweight + &
                       D0* (-1.D0+pc_p_loc_p(jng))+gdz* (1.D0- &
                       upweight)*d_p_loc_p(jng)*avgmw_loc_p(jng))
      !    print *,'liquid devq::',devq(1,2),kvr_p_loc_p(jng),pc_loc_p(jng),&
      !             pc_p_loc_p(jng),d_p_loc_p(jng),ng,jng
          devq(2,2)=Dq*(kvr_t_loc_p(jng)*dphi*mnweight + D0* &
                    (pc_t_loc_p(jng)+ gdz*(1.D0-upweight)*   &
                    avgmw_loc_p(jng) *  d_t_loc_p(jng)))
      !    print *,'liquid devq::',devq(2,2)
          devq(3,2)=Dq*(kvr_c_loc_p(jng)*dphi*mnweight &
                     +D0*(pc_c_loc_p(jng)+ gdz*(1.D0-upweight)* &
                     (avgmw_loc_p(jng) *d_c_loc_p(jng)+ &
                     ddensity_loc_p(jng)*avgmw_c_loc_p(jng))))
      !    print *,'liquid devq::',devq(3,2)
          devq(4,2)=Dq*(kvr_s_loc_p(jng)*dphi*mnweight + &
                     D0*pc_s_loc_p(jng))
      !    print *,'liquid devq::',devq(4,2)
          devq=devq* grid%areabc(nc)
      !    print *,'liquid devq::',devq
! pressure equation:: Boundary do not have equation so ra(1,1:4)=0
        ra(1,5)=ra(1,5) + devq(1,2)*density_ave + q*(1.D0-upweight)*d_p_loc_p(jng)
        ra(1,6)=ra(1,6) + devq(2,2)*density_ave + q*(1.D0-upweight)*d_t_loc_p(jng)
        ra(1,7)=ra(1,7) + devq(3,2)*density_ave + q*(1.D0-upweight)*d_c_loc_p(jng)
        ra(1,8)=ra(1,8) + devq(4,2)*density_ave + q*(1.D0-upweight)*d_s_loc_p(jng)

! energy equation
          ra(2,5)= ra(2,5) + uhh*density_ave*devq(1,2) &
                   + q * uhh*d_p_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_p_loc_p(jng)* mnweight
       
           ra(2,6)= ra(2,6)+ uhh*density_ave*devq(2,2) &
                   + q * uhh*d_t_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_t_loc_p(jng)* mnweight
          
           ra(2,7)= ra(2,7)+uhh*density_ave*devq(3,2) &
                   + q * uhh*d_c_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_c_loc_p(jng)* mnweight
       
           ra(2,8)= ra(2,8)+ uhh*density_ave*devq(4,2) &
                   + q * uhh*d_s_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_s_loc_p(jng)* mnweight

! Concentration equation
           tempvar=0.D0
     !Contribution from convection    

           ibc_hencoeff=2+(j-1)*grid%nspec
           i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
           i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(ng-1)*grid%nphase*&
                grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp

           cw1 =1.D0-(grid%concbc(nc))*grid%hen_bc(ibc_hencoeff)
           cw2 =1.D0-(CCONC_LOC(ng))*Hen_loc_p(i2_hencoeff)

           if(mu==0)then
              cwu=cw1
           else
              cwu=cw2  
           end if

           cw2_p= -(CCONC_LOC(ng))*Hen_p_loc_p(i2_hencoeff)
           cw2_t= -(CCONC_LOC(ng))*Hen_t_loc_p(i2_hencoeff)
           cw2_c= -Hen_loc_p(i2_hencoeff)-(CCONC_LOC(ng))*Hen_c_loc_p(i2_hencoeff_dc)
           cw2_s= -(CCONC_LOC(ng))*Hen_s_loc_p(i2_hencoeff)

           tempvar(5) = d_p_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(1,2) * cwu   &
                + density_ave * q * mnweight *cw2_p

           tempvar(6) =  d_t_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(2,2) * cwu   &
                + density_ave * q * mnweight * cw2_t
    
           tempvar(7) = d_c_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(3,2) * cwu   &
                + density_ave * q * mnweight*cw2_c
           tempvar(8) = d_s_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(4,2) * cwu   &
                + density_ave * q * mnweight * cw2_s

           ra(4,5:8)=ra(4,5:8)+tempvar(5:8)
    !       ra(5,5:8)=ra(5,5:8)+tempvar(5:8)
        end if
     end if
!print *,'     In 2ph Jacobian ::  BC 1 : liquid conv '
 
 !print *,ra   
! ************** For Gas phase Flux ********************
     if((grid%sgbc(nc)> eps).or.(SSATG_LOC(ng)>eps))then
    
            j=2
            jng = j + (ng-1) * grid%nphase 

            upweight=1.d0
            if(grid%sgbc(nc)<eps)then
              upweight=0.d0
            endif
            density_ave = upweight* grid%density_bc(j)+ &
                        (1.D0-upweight)*ddensity_loc_p(jng)  

            gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                 (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                 * grid%gravity * grid%delzbc(nc)
            dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc)- gravity)
            
          if(dphi>=0.D0)then
             mu=0 ! BC
             ukvr=grid%kvr_bc(j)
             uhh=grid%hh_bc(j)
             uconc=grid%concbc(nc)
            else 
             mu=ng
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc=CCONC_LOC(ng)
          endif

           if((ukvr*Dq)>floweps)then 

              D0=ukvr 

              v_darcy =  Dq * D0 * dphi
              q = v_darcy * grid%areabc(nc)
              !q = 0.d0
              gdz=grid%gravity * grid%delzbc(nc)
                              
              if(mu==0)then   
                 mbweight = 1.D0
                 mnweight = 0.D0
              else
                 mbweight = 0.D0
                 mnweight = 1.D0
              endif

! now 1 refer to boundary
!     2          node contains BC
! So do not take derivative to index 1
          devq=0.D0    
          devq(1,2) = Dq*(kvr_p_loc_p(jng)*dphi*mnweight + &
                       D0* (-1.D0+ gdz* (1.D0- &
                       upweight)*d_p_loc_p(jng)*avgmw_loc_p(jng)))

          devq(2,2)=Dq*(kvr_t_loc_p(jng)*dphi*mnweight + D0* &
                    ( gdz*(1.D0-upweight)*   &
                    avgmw_loc_p(jng) *  d_t_loc_p(jng)))

          devq(3,2)=Dq*(kvr_c_loc_p(jng)*dphi*mnweight &
                     +D0*(pc_c_loc_p(jng)+ gdz*(1.D0-upweight)* &
                     (avgmw_loc_p(jng) *d_c_loc_p(jng)+ &
                     ddensity_loc_p(jng)*avgmw_c_loc_p(jng))))

          devq(4,2)=Dq*(kvr_s_loc_p(jng)*dphi*mnweight + &
                     D0*pc_s_loc_p(jng))

          devq=devq* grid%areabc(nc)
         

! pressure equation:: Boundary do not have equation so ra(1,1:4)=0
          ra(1,5)=ra(1,5) + devq(1,2)*density_ave+q*(1.D0-upweight)*d_p_loc_p(jng)
          ra(1,6)=ra(1,6) + devq(2,2)*density_ave+q*(1.D0-upweight)*d_t_loc_p(jng)
          ra(1,7)=ra(1,7) + devq(3,2)*density_ave+q*(1.D0-upweight)*d_c_loc_p(jng)
          ra(1,8)=ra(1,8) + devq(4,2)*density_ave+q*(1.D0-upweight)*d_s_loc_p(jng)

! energy equation
           ra(2,5)= ra(2,5) + uhh*density_ave*devq(1,2) &
                   + q * uhh*d_p_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_p_loc_p(jng)* mnweight
       
           ra(2,6)= ra(2,6)+ uhh*density_ave*devq(2,2) &
                   + q * uhh*d_t_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_t_loc_p(jng)* mnweight
          
           ra(2,7)= ra(2,7)+uhh*density_ave*devq(3,2) &
                   + q * uhh*d_c_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_c_loc_p(jng)* mnweight
       
           ra(2,8)= ra(2,8)+ uhh*density_ave*devq(4,2) &
                   + q * uhh*d_s_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_s_loc_p(jng)* mnweight


           uconc=1.D0-uconc
           tempvar(5)= uconc*density_ave*devq(1,2) &
                   + q *uconc *d_p_loc_p(jng)*(1.D0-upweight)
                  
           tempvar(6)= uconc*density_ave*devq(2,2) &
                    + q * uconc *d_t_loc_p(jng)*(1.D0-upweight)
           tempvar(7)= uconc*density_ave*devq(3,2) &
                   + q * uconc *d_c_loc_p(jng)*(1.D0-upweight) &
                   - q * density_ave * mnweight    
           tempvar(8)= uconc*density_ave*devq(4,2) &
                   + q * uconc *d_s_loc_p(jng)*(1.D0-upweight)
     
  
           ra(3,5:8)=ra(3,5:8)+tempvar(5:8)
! Saturation equation
           ra(4,5:8)=ra(4,5:8)+tempvar(5:8) 
!print *,'     In 2ph Jacobian ::  BC 1 : gas conv '
        end if
     end if
!     print *,'BC xmol,s ',nc,uconc
!print*,ra

  case(2)  
     if((dabs(grid%velocitybc(1,nc))+dabs(grid%velocitybc(2,nc)))>floweps)then
     flux =0.; fluxh=0.; fluxv=0.;
     tempvar=0.D0

     do j=1,grid%nphase
        jng = j + (ng-1) * grid%nphase
        v_darcy = grid%velocitybc(j,nc)
        if(v_darcy < 0.)then 
           q = v_darcy * ddensity_loc_p(jng) * grid%areabc(nc)
           !q = 0.d0
           !v_darcy = 0.d0
! pressure equation
           ra(1,5)=ra(1,5)+ v_darcy *grid%areabc(nc)*d_p_loc_p(jng)
           ra(1,6)=ra(1,6)+ v_darcy *grid%areabc(nc)*d_t_loc_p(jng)
           ra(1,7)=ra(1,7)+ v_darcy *grid%areabc(nc)*d_c_loc_p(jng)
           ra(1,8)=ra(1,8)+ v_darcy *grid%areabc(nc)*d_s_loc_p(jng)
! energy equation
           ra(2,5)=ra(2,5)+ v_darcy *grid%areabc(nc)*(d_p_loc_p(jng)* &
                    hh_loc_p(jng)+ddensity_loc_p(jng)*h_p_loc_p(jng))
           ra(2,6)=ra(2,6)+ v_darcy *grid%areabc(nc)*(d_t_loc_p(jng)* &
                    hh_loc_p(jng)+ddensity_loc_p(jng)*h_t_loc_p(jng))
           ra(2,7)=ra(2,7)+ v_darcy *grid%areabc(nc)*(d_c_loc_p(jng)* &
                    hh_loc_p(jng)+ddensity_loc_p(jng)*h_c_loc_p(jng))
           ra(2,8)=ra(2,8)+ v_darcy *grid%areabc(nc)*(d_s_loc_p(jng)* &
                    hh_loc_p(jng)+ddensity_loc_p(jng)*h_s_loc_p(jng))

!Concentration equation
           uconc=1.D0-CCONC_LOC(ng)
           select case(j)
           case(2)
              tempvar(5)=tempvar(5) + v_darcy *grid%areabc(nc)*(d_p_loc_p(jng)* uconc)
              tempvar(6)=tempvar(6) + v_darcy *grid%areabc(nc)*(d_t_loc_p(jng)* uconc)
              tempvar(7)= tempvar(7)+ v_darcy *grid%areabc(nc)*(d_c_loc_p(jng)* &
                   uconc - ddensity_loc_p(jng))
              tempvar(6)=tempvar(6) + v_darcy *grid%areabc(nc)*(d_s_loc_p(jng)* uconc)

            ra(3,5:8)=ra(3,5:8) + tempvar(5:8)
            ra(3,5:8)=ra(3,5:8) + tempvar(5:8)
           case(1)
              i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
              i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(ng-1)*grid%nphase*&
                  grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp

              cw =1.D0-(CCONC_LOC(ng))* hen_loc_p(i2_hencoeff)
              cw_p=-(CCONC_LOC(ng))* hen_p_loc_p(i2_hencoeff)
              cw_t=-(CCONC_LOC(ng))* hen_t_loc_p(i2_hencoeff)
              cw_c=-hen_loc_p(ng)-(CCONC_LOC(ng))*hen_c_loc_p(i2_hencoeff_dc)
              cw_s=-(CCONC_LOC(ng))* hen_s_loc_p(i2_hencoeff)

              tempvar(5)=v_darcy *grid%areabc(nc)*(d_p_loc_p(jng)* &
                      cw + ddensity_loc_p(jng)*cw_p)
              tempvar(6)=v_darcy *grid%areabc(nc)*(d_t_loc_p(jng)* &
                   cw + ddensity_loc_p(jng)*cw_t)
              tempvar(7)=v_darcy *grid%areabc(nc)*(d_c_loc_p(jng)* &
                   cw + ddensity_loc_p(jng)*cw_c)
              tempvar(8)=v_darcy *grid%areabc(nc)*(d_s_loc_p(jng)* &
                   cw + ddensity_loc_p(jng)*cw_s)
              
              ra(4,5:8)=ra(4,5:8)+tempvar(5:8)
              ra(5,5:8)=ra(4,5:8)+tempvar(5:8)

           end select
         end if
      end do
    endif
   case(3)   
  
     Dq= perm1 / grid%distbc(nc)

      ! Liquid phase   
         if((1.D0-SSATG_LOC(ng))>grid%swir(int(icap_loc_p(ng))))then
          
            j=1
            jng = j + (ng-1) * grid%nphase 

         
            upweight=0.D0
          
            density_ave = upweight* grid%density_bc(j)+ &
                        (1.D0-upweight)*ddensity_loc_p(jng)  

            gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                 (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                 * grid%gravity * grid%delzbc(nc)
            dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc) & 
                 - pc_loc_p(jng) + grid%pc_bc(j)   &
                 - gravity)
            
       
             mu=ng
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc=CCONC_LOC(ng)
             

           if((ukvr*Dq)>floweps)then 

              D0=ukvr 

              v_darcy =  Dq * D0 * dphi
              q = v_darcy * grid%areabc(nc)
              !q = 0.d0
              gdz=grid%gravity * grid%delzbc(nc)
                              
      
                 mbweight = 0.D0
                 mnweight = 1.D0
      
! now 1 refer to boundary
!     2          node contains BC
! So do not take derivative to index 1
          devq=0.D0    
          devq(1,2) = Dq*(kvr_p_loc_p(jng)*dphi*mnweight + &
                       D0* (-1.D0+pc_p_loc_p(jng))+gdz* (1.D0- &
                       upweight)*d_p_loc_p(jng)*avgmw_loc_p(jng))

          devq(2,2)=Dq*(kvr_t_loc_p(jng)*dphi*mnweight + D0* &
                    (pc_t_loc_p(jng)+ gdz*(1.D0-upweight)*   &
                    avgmw_loc_p(jng) *  d_t_loc_p(jng)))

          devq(3,2)=Dq*(kvr_c_loc_p(jng)*dphi*mnweight &
                     +D0*(pc_c_loc_p(jng)+ gdz*(1.D0-upweight)* &
                     (avgmw_loc_p(jng) *d_c_loc_p(jng)+ &
                     ddensity_loc_p(jng)*avgmw_c_loc_p(jng))))

          devq(4,2)=Dq*(kvr_s_loc_p(jng)*dphi*mnweight + &
                     D0*pc_s_loc_p(jng))
          devq=devq* grid%areabc(nc)
         
! pressure equation:: Boundary do not have equation so ra(1,1:4)=0
          ra(1,5)=ra(1,5) + devq(1,2)*density_ave + q*(1.D0-upweight)*d_p_loc_p(jng)
          ra(1,6)=ra(1,6) + devq(2,2)*density_ave + q*(1.D0-upweight)*d_t_loc_p(jng)
          ra(1,7)=ra(1,7) + devq(3,2)*density_ave + q*(1.D0-upweight)*d_c_loc_p(jng)
          ra(1,8)=ra(1,8) + devq(4,2)*density_ave + q*(1.D0-upweight)*d_s_loc_p(jng)

! energy equation
          ra(2,5)= ra(2,5) + uhh*density_ave*devq(1,2) &
                   + q * uhh*d_p_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_p_loc_p(jng)* mnweight
       
           ra(2,6)= ra(2,6)+ uhh*density_ave*devq(2,2) &
                   + q * uhh*d_t_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_t_loc_p(jng)* mnweight
          
           ra(2,7)= ra(2,7)+uhh*density_ave*devq(3,2) &
                   + q * uhh*d_c_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_c_loc_p(jng)* mnweight
       
           ra(2,8)= ra(2,8)+ uhh*density_ave*devq(4,2) &
                   + q * uhh*d_s_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_s_loc_p(jng)* mnweight

! Concentration equation
!           ibc_hencoeff=1+(j-1)*grid%nspec
           i2_hencoeff=2+(j-1)*grid%nspec+(ng-1)*grid%nphase*grid%nspec
           i2_hencoeff_dc=2+(j-1)*grid%npricomp*grid%nspec+(ng-1)*grid%nphase*&
                grid%nspec*grid%npricomp ! +(jc-1)*grid%npricomp


           cw2 =1.D0-(CCONC_LOC(ng))*Hen_loc_p(i2_hencoeff)
           cw2_p= -(CCONC_LOC(ng))*Hen_p_loc_p(i2_hencoeff)
           cw2_t= -(CCONC_LOC(ng))*Hen_t_loc_p(i2_hencoeff)
           cw2_c= -Hen_loc_p(i2_hencoeff)-(1.D0-CCONC_LOC(ng))*Hen_c_loc_p(i2_hencoeff_dc)
           cw2_s= -(CCONC_LOC(ng))*Hen_s_loc_p(i2_hencoeff)

           tempvar=0.D0
         tempvar(5)= d_p_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(1,2) * cwu   &
                + density_ave * q * mnweight *cw2_p

         tempvar(6)= d_t_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(2,2) * cwu   &
                + density_ave * q * mnweight * cw2_t
    
         tempvar(7)= d_c_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(3,2) * cwu   &
                + density_ave * q * mnweight*cw2_c
         tempvar(8)= d_s_loc_p(jng)*(1.D0-upweight) *q * cwu &
                + density_ave *devq(4,2) * cwu   &
                + density_ave * q * mnweight * cw2_s

          ra(4,5:8)= ra(4,5:8)+ tempvar(5:8)
    !      ra(5,5:8)= ra(5,5:8)+ tempvar(5:8)

        end if
     end if

! ************** For Gas phase Flux ********************
     if(SSATG_LOC(ng)>eps)then
    
            j=2
            jng = j + (ng-1) * grid%nphase 

               upweight=0.D0
            
            density_ave =ddensity_loc_p(jng)  

            gravity = (upweight * grid%density_bc(j)*grid%avgmw_bc(j) + &
                 (1.D0-upweight)*ddensity_loc_p(jng)*avgmw_loc_p(jng)) &
                 * grid%gravity * grid%delzbc(nc)
            dphi= -(PPRESSURE_LOC(ng) - grid%pressurebc(2,nc)- gravity)
       
             mu=ng
             ukvr=kvr_loc_p(jng)
             uhh=hh_loc_p(jng)
             uconc= CCONC_LOC(ng)

       
           if((ukvr*Dq)>floweps)then 

              D0=ukvr 

              v_darcy =  Dq * D0 * dphi
              q = v_darcy * grid%areabc(nc)
              !q = 0.d0
              gdz=grid%gravity * grid%delzbc(nc)
          
                 mbweight = 0.D0
                 mnweight = 1.D0
             
! now 1 refer to boundary
!     2          node contains BC
! So do not take derivative to index 1
          devq=0.D0    
          devq(1,2) = Dq*(kvr_p_loc_p(jng)*dphi*mnweight + &
                       D0* (-1.+ gdz* (1.D0- &
                       upweight)*d_p_loc_p(jng)*avgmw_loc_p(jng)))

          devq(2,2)=Dq*(kvr_t_loc_p(jng)*dphi*mnweight + D0* &
                    ( gdz*(1.D0-upweight)*   &
                    avgmw_loc_p(jng) *  d_t_loc_p(jng)))

          devq(3,2)=Dq*(kvr_c_loc_p(jng)*dphi*mnweight &
                     +D0*(pc_c_loc_p(jng)+ gdz*(1.D0-upweight)* &
                     (avgmw_loc_p(jng) *d_c_loc_p(jng)+ &
                     ddensity_loc_p(jng)*avgmw_c_loc_p(jng))))

          devq(4,2)=Dq*(kvr_s_loc_p(jng)*dphi*mnweight + &
                     D0*pc_s_loc_p(jng))

          devq=devq* grid%areabc(nc)

! pressure equation:: Boundary do not have equation so ra(1,1:4)=0
          ra(1,5)=ra(1,5) + devq(1,2)*density_ave+q*(1.D0-upweight)*d_p_loc_p(jng)
          ra(1,6)=ra(1,6) + devq(2,2)*density_ave+q*(1.D0-upweight)*d_t_loc_p(jng)
          ra(1,7)=ra(1,7) + devq(3,2)*density_ave+q*(1.D0-upweight)*d_c_loc_p(jng)
          ra(1,8)=ra(1,8) + devq(4,2)*density_ave+q*(1.D0-upweight)*d_s_loc_p(jng)

! energy equation
           ra(2,5)= ra(2,5) + uhh*density_ave*devq(1,2) &
                   + q * uhh*d_p_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_p_loc_p(jng)* mnweight
       
           ra(2,6)= ra(2,6)+ uhh*density_ave*devq(2,2) &
                   + q * uhh*d_t_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_t_loc_p(jng)* mnweight
          
           ra(2,7)= ra(2,7)+uhh*density_ave*devq(3,2) &
                   + q * uhh*d_c_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_c_loc_p(jng)* mnweight
       
           ra(2,8)= ra(2,8)+ uhh*density_ave*devq(4,2) &
                   + q * uhh*d_s_loc_p(jng)*(1.D0-upweight) &
                   + q * density_ave * h_s_loc_p(jng)* mnweight



! Saturation equation
        tempvar=0.D0
     !Contribution from convection    
           uconc=1.D0-uconc
           tempvar(5)= uconc*density_ave*devq(1,2) &
                   + q *uconc *d_p_loc_p(jng)*(1.D0-upweight)
                  
           tempvar(6)= uconc*density_ave*devq(2,2) &
                    + q * uconc *d_t_loc_p(jng)*(1.D0-upweight)
           tempvar(7)= uconc*density_ave*devq(3,2) &
                   + q * uconc *d_c_loc_p(jng)*(1.D0-upweight) &
                   - q * density_ave * mnweight    
           tempvar(8)= uconc*density_ave*devq(4,2) &
                   + q * uconc *d_s_loc_p(jng)*(1.D0-upweight)
     
  


           ra(3,5:8)=ra(3,5:8)+tempvar(5:8)
! Saturation equation
           ra(4,5:8)=ra(4,5:8)+tempvar(5:8)


  end if
end if

       
end select

!    if(SSATG_LOC(ng)<eps) ra(:,8)=0.D0

tempvar=0.D0
 ! if(SSATG_LOC(ng)<=eps)then
 !   ra(3,4:8)= ra(5,4:8)
! endif
    tempvar(4:8)=ra(3,4:8)
     ra(3,4:8)=ra(4,4:8)
     ra(4,4:8)=tempvar(4:8)
!  end if

   blkmat11 = 0.D0;
   do ii=0,3
        do jj=4,7
           select case(ii)
              case(3)
                 if(iiphas==6)then
           if (grid%iblkfmt == 0) then
             call MatSetValuesLocal(A,1,p1+ii,1,p1+jj-4,-ra(ii+1,jj+1),ADD_VALUES,ierr)
!                    else
 !                      call MatSetValuesLocal(A,1,p1+ii,1,p1+jj-4,-ra(5,jj+1),ADD_VALUES,ierr)
                    else
              blkmat11(ii+1,jj-4+1) =blkmat11(ii+1,jj-4+1) - ra(ii+1,jj+1)
            endif
           end if
              case default
          if (grid%iblkfmt == 0) then
                   call MatSetValuesLocal(A,1,p1+ii,1,p1+jj-4,-ra(ii+1,jj+1),ADD_VALUES,ierr)
         else
           blkmat11(ii+1,jj-4+1) =blkmat11(ii+1,jj-4+1) - ra(ii+1,jj+1)
         endif  
              end select
        enddo
     enddo
   
  if (grid%iblkfmt == 1) then
    call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat11,ADD_VALUES,ierr)
  endif  
   
   
 
  end do


  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)

! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%avgmw_loc, avgmw_loc_p, ierr)
!  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(grid%uu, uu_p, ierr)
  call VecRestoreArrayF90(grid%df_loc, df_loc_p, ierr)
  call VecRestoreArrayF90(grid%hen_loc, hen_loc_p, ierr)
  call VecRestoreArrayF90(grid%pcw_loc, pc_loc_p, ierr)
  call VecRestoreArrayF90(grid%kvr_loc, kvr_loc_p, ierr)


    call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_c_loc, d_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_s_loc, d_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%avgmw_c_loc,avgmw_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_c_loc, h_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_s_loc, h_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%u_p, u_p_p, ierr)
    call VecRestoreArrayF90(grid%u_t, u_t_p, ierr)
    call VecRestoreArrayF90(grid%u_c, u_c_p, ierr)
    call VecRestoreArrayF90(grid%u_s, u_s_p, ierr)
    call VecRestoreArrayF90(grid%df_p_loc, df_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_t_loc, df_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_c_loc, df_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%df_s_loc, df_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_p_loc, hen_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_t_loc, hen_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_c_loc, hen_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%hen_s_loc, hen_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_p_loc, pc_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_t_loc, pc_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_c_loc, pc_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%pc_s_loc, pc_s_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_p_loc, kvr_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_t_loc, kvr_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_c_loc, kvr_c_loc_p, ierr)
    call VecRestoreArrayF90(grid%kvr_s_loc, kvr_s_loc_p, ierr)
   
    call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
    call VecRestoreArrayF90(grid%icap_loc, icap_loc_p, ierr)
    call VecRestoreArrayF90(grid%iphas_loc, iphase_loc_p, ierr)

if (grid%rk > 0.d0) then
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif

  

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

!  B = A
 
 !call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB, ierr)
   

! call MatView(A, PETSC_VIEWER_STDOUT_WORLD,ierr)
! stop

end subroutine TTPHASEJacobin



 subroutine pflow_2phase_initaccum(grid)
 
  use mixture_module  
   implicit none
    type(pflowGrid) :: grid 

 
  integer :: ierr
  integer :: n
  integer :: i, j, jn
  integer :: p1, t1, c1, s1
  integer :: ii1,ii2, iicap

  PetscScalar, pointer ::accum_p(:),yy_p(:)
  
  PetscScalar, pointer ::  porosity_p(:), volume_p(:), &
                           density_p(:), avgmw_p(:), u_p(:), &
                           h_p(:), df_p(:), hen_p(:),&
                           pc_p(:), kvr_p(:),&
                           ithrm_p(:),icap_p(:),iphase_p(:) 

!  integer, pointer ::iphase_p(:)
  
  real*8 :: sat_pressure, pvol, satw  ! Saturation pressure of water.
 
  real*8 :: xxlw,xxla,xxga,xxgw,eengl,eengg,acc

  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecGetArrayF90(grid%h, h_p, ierr)
  call VecGetArrayF90(grid%u, u_p, ierr)
  call VecGetArrayF90(grid%pcw, pc_p, ierr)
  call VecGetArrayF90(grid%kvr, kvr_p, ierr)
  call VecGetArrayF90(grid%hen, hen_p,ierr) 
  call VecGetArrayF90(grid%df, df_p,ierr) 
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)

 do n = 1, grid%nlmax
        jn = 1 + (n-1)*grid%nphase
        ii1=1+(n-1)*grid%nphase; ii2=n*grid%nphase
        iicap=int(icap_p(n))

      call mixture_eos_noderiv (PRESSURE(n),TEMP(n),CONC(n),SATG(n),&
        grid%scale,grid%nphase,grid%nspec,grid%npricomp,&
        iicap,grid%swir(iicap),grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),&
        grid%pcwmax(iicap),grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
        density_p(ii1:ii2),avgmw_p(ii1:ii2),& 
        h_p(ii1:ii2),u_p(ii1:ii2),&
        df_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        grid%xxphi_co2(n),&
        pc_p(ii1:ii2),kvr_p(ii1:ii2),&
        ierr,grid%itable)

!     if(n==1)then
!     if(n==3585)then
!       print *,'init acc: ',n,ii1,jn,iicap,PRESSURE(n),&
!                  TEMP(n),CONC(n),SATG(n),density_p(ii1:ii2)

         !print *,  grid%scale,grid%nphase,grid%nspec,grid%npricomp
         !print *,'Cap  ',  iicap, grid%swir(iicap), grid%lambda(iicap),&
!                       grid%alpha(iicap),grid%pckrm(iicap),&
!                       grid%pcwmax(iicap)
        !print *, 'Ps  ', sat_pressure
        !print *, 'Pressure', PRESSURE(n)
        !print *, 'den   ', density_p(ii1:ii2),avgmw_p(ii1:ii2)
        !print *, 'H  U    ',   h_p(ii1:ii2),u_p(ii1:ii2)
        !print *,'diff   ', df_p(1+(n-1)*grid%nphase*grid%nspec:&
!                            n*grid%nphase*grid%nspec)
        !print *,'Henry  ',hen_p(1+(n-1)*grid%nphase*grid%nspec:&
!                            n*grid%nphase*grid%nspec)
        !print *,'pc  ', pc_p(ii1:ii2)
        !print *,'kvr   ', kvr_p(ii1:ii2)
!     endif

 enddo

 
!---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
  !  ng = grid%nL2G(n)   ! corresponding ghost index
    p1 = 1 + (n-1)*grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1
    s1 = c1 + grid%npricomp

    pvol=volume_p(n)*porosity_p(n)
    SATW=1.D0-SATG(n)

! Pressure equation accumulation term
    acc = 0.d0
!    do j = 1, grid%nphase
      j=1
      jn = j + (n-1)*grid%nphase
  
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
      acc = acc +   pvol*(density_p(jn)*SATW     &
              +density_p(1+jn)*SATG(n))   
                           
!    enddo
    accum_p(p1) = acc
    !print *, 'P accum ', n,  accum_p(p1) 

! Heat equation accumulation term
    i = ithrm_p(n)
!    j = grid%jh2o
    j=1
    jn = j+(n-1)*grid%nphase
 
    eengl = density_p(jn)*u_p(jn)
    eengg = density_p(1+jn)*u_p(1+jn)
    eengl=eengl * SATW
    eengg=eengg * SATG(n)
 
 !   engl = density_p(jn) * h_p(jn) - grid%scale*PRESSURE(n)
 !   engg = density_p(1+jn) * h_p(1+jn) - grid%scale*PRESSURE(n)  
 !   engl=engl * ( 1.D0 - SATG(n) )
 !   engg=engg * SATG(n)

    accum_p(t1) = pvol * (eengl+eengg) &     ! fluid
                + (1.d0-porosity_p(n)) * grid%dencpr(i) * &
                TEMP(n)*volume_p(n)   ! Rock 
               
    !print *, 'T accum ', n,  accum_p(t1) 
   
! Concentration equation accumulation term
! remember we had choosed the X(water) in liquid as primary variable
 
      !print *, 'C accum ', n,  accum_p(c1)          

! Saturation equation accumulation term
   ! call Henry_coeff(PPRESSURE_LOC(ng),TTEMP_LOC(ng),henrycoeff)
    xxga = CONC(n)
    xxgw = 1.D0-xxga
    xxla =  xxga*hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec)
    xxlw = 1.D0 - xxla
!    yylw = CCONC(n)
!    yyla = 1.D0-yylw
!    yyga = henry_p(n) * yyla
!    yygw = 1.D0 - yyga
!    print *,'initaccum', xxla,hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec)
  
  accum_p(c1) = pvol * density_p(jn+1) * xxgw * SATG(n) &
               * grid%ret
  accum_p(s1) = pvol * (density_p(jn) * xxlw * SATW) * grid%ret + accum_p(c1)
! print *, 'init S accum ', n,  accum_p(s1) 

! print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec),kvr_p(jn),kvr_p(jn+1)

 end do

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%accum, accum_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecRestoreArrayF90(grid%h, h_p, ierr)
  call VecRestoreArrayF90(grid%u, u_p, ierr)
  call VecRestoreArrayF90(grid%pcw, pc_p, ierr)
  call VecRestoreArrayF90(grid%kvr, kvr_p, ierr)
  call VecRestoreArrayF90(grid%hen, hen_p,ierr) 
  call VecRestoreArrayF90(grid%df, df_p,ierr)
  call VecRestoreArrayF90(grid%icap, icap_p, ierr) 
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  
end subroutine pflow_2phase_initaccum

  subroutine pflow_update_2phase(grid)
    use mixture_module  
    use water_eos_module
    implicit none
    type(pflowGrid) :: grid 

  
    PetscScalar, pointer :: t_p(:),p_p(:),c_p(:),s_p(:),cc_p(:)
    PetscScalar, pointer :: yy_p(:),pc_p(:),hen_p(:)
                        
    integer :: n, ii1, ii2, jn
    integer :: ierr
 
    real*8 :: xxlw,xxla,xxga,xxgw

    if(grid%nphase>1) call pflow_2phase_massbal(grid)
  ! if (grid%rk > 0.d0) call Rock_Change(grid)
    call  TTPhase_Update(grid%xx,grid)
!   call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
!   call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
!   call VecGetArrayF90(grid%iphas, iphase_p, ierr); CHKERRQ(ierr)
! do n = 1, grid%nlmax
!   jn = (n-1)*grid%ndof

   ! if( (xx_p(jn+4)<satcuteps) .and. (xx_p(jn+4)<yy_p(jn+4)) ) xx_p(jn+4)=0.D0
   ! if( (xx_p(jn+4)>1.D0-satcuteps) .and. (xx_p(jn+4)>yy_p(jn+4)))&
   !      xx_p(jn+4)=1.D0 

!   if(  iphase_p(n)==2)then 
!      iphase_p(n) = 2
!   elseif (iphase_p(n) == 4)then
!      iphase_p(n) = 4
!   else
!      iphase_p(n)=6
!    endif

 ! if(xx_p(jn+4)<eps  iphase_p(n) = 2
 !   if(xx_p(jn+4)>(1.D0-eps)) iphase_p(n) = 4


!    if( (xx_p(jn+4)>eps) .and. ( yy_p(jn+4)<=eps) &
!         .and. (xx_p(jn+4)>yy_p(jn+4)))then
!       xx_p(jn+3)=1.D0
!       print *,'init break  ',n
!    endif
 
   ! if( (xx_p(jn+4)<=1.D0-satcuteps) .and. (yy_p(jn+4)>=1.D0-satcuteps))then
   !       xx_p(jn+3)=1.D0-sat_pressure/xx_p(jn+1)
   !       print *,'update adjust L->G',n,xx_p(jn+1:jn+4)
   !       endif

!print *,'updating ',n,xx_p(jn+1:jn+4),sat_pressure,1.D0-sat_pressure/xx_p(jn+1)
!print *,'         ',n, yy_p(jn+1:jn+4),1.D0-satcuteps

!  call PSAT(xx_p(jn+2), sat_pressure, ierr)
 !   if((xx_p(jn+4)<satcuteps) .and. xx_p(jn+1)<sat_pressure)then
 !      xx_p(jn+4)=satcuteps*1.001D0
 !      xx_p(jn+3)=0.D0
 !      print *,'update adjust L->G',n,xx_p(jn+4), xx_p(jn+1)
 !   end if
 !   if((xx_p(jn+4)>1.D0-satcuteps) .and.(xx_p(jn+1)*(1.D0-xx_p(jn+3)))>sat_pressure)then
 !      xx_p(jn+4)=(1.D0-satcuteps)* 0.999D07
 !      xx_p(jn+3)=1.D0-sat_pressure/xx_p(jn+1)
 !      print *,'update adjust G->L',n,xx_p(jn+4), xx_p(jn+1)
 !   end if

  !  if( (xx_p(jn+3)>1.D0-eps).and.(xx_p(jn+3)<1.D0+eps))xx_p(jn+3)=1.D0-eps
  !  if( (xx_p(jn+3)>-eps).and.(xx_p(jn+3)<eps))xx_p(jn+3)=0.D0
! if(n==grid%nlmax/2)print*,n,xx_p(jn+1:jn+4)
!enddo


! call VecRestoreArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
!  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
!  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)

  call VecCopy(grid%xx, grid%yy, ierr)   
  call VecCopy(grid%iphas, grid%iphas_old, ierr)   
 
  call  pflow_2phase_initaccum(grid)

  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%pressure, p_p, ierr)
  call VecGetArrayF90(grid%temp, t_p, ierr)
  call VecGetArrayF90(grid%xmol, c_p, ierr)
  call VecGetArrayF90(grid%sat, s_p, ierr)
  call VecGetArrayF90(grid%pcw, pc_p, ierr)
  call VecGetArrayF90(grid%hen, hen_p, ierr)
  call VecGetArrayF90(grid%conc,cc_p, ierr)



 do n = 1, grid%nlmax
   
   jn = 1 + (n-1)*grid%nphase
        ii1 = jn
        ii2=jn + grid%nphase -1


      p_p(jn)  = PRESSURE(n)- pc_p(ii1)
      p_p(jn+1)= PRESSURE(n)
      t_p(n) = TEMP(n)
      c_p(jn+1) = CONC(n)
      xxga = CONC(n)
      xxgw = 1.D0-xxga
      xxla =  xxga*hen_p(2+(n-1)*grid%nphase*grid%nspec)
      xxlw = 1.D0 - xxla
      c_p(jn) = xxla 
      cc_p(n)=xxga
      s_p(jn+1)= SATG(n)
      s_p(jn)= 1.D0 - SATG(n)
  !    print *, 'updating...'
  !    print *, PRESSURE(n), pc_p(ii1), TEMP(n),CONC(n), SATG(n)
   enddo
!   print *,' 2 ph update : hen', hen_p(1+(grid%nlmax/2-1)*grid%nphase*grid%nspec)

 call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%pressure, p_p, ierr)
  call VecRestoreArrayF90(grid%temp, t_p, ierr)
  call VecRestoreArrayF90(grid%xmol, c_p, ierr)
  call VecRestoreArrayF90(grid%sat, s_p, ierr)
 call VecRestoreArrayF90(grid%pcw, pc_p, ierr)
  call VecRestoreArrayF90(grid%hen, hen_p, ierr)
  call VecRestoreArrayF90(grid%conc, cc_p,ierr)

 end subroutine pflow_update_2phase






  subroutine pflow_2phase_initadj(grid)
 
! running this subroutine will override the xmol data for initial condition in pflow.in 

  use mixture_module  
  implicit none
  type(pflowGrid) :: grid 

 
  integer :: ierr
  integer :: n, nc
  integer :: ibc,jn
  integer :: m
  integer :: ii1,ii2,iicap
  integer :: i1_hencoeff
 

  PetscScalar, pointer ::  p_p(:), t_p(:), x_p(:), s_p(:),iphase_p(:)
  PetscScalar, pointer ::  porosity_p(:), volume_p(:), &
                           density_p(:), avgmw_p(:), u_p(:), &
                           h_p(:), df_p(:), hen_p(:),&
                           pc_p(:), kvr_p(:),&
                           ithrm_p(:),icap_p(:) 

!  integer, pointer ::iphase_p(:)
  
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: xlw,xla,xga,xgw
  
! real*8 :: temp1
  real*8, parameter :: Rg=8.31415D0

  call VecGetArrayF90(grid%pressure, p_p, ierr)
  call VecGetArrayF90(grid%temp, t_p, ierr)
  call VecGetArrayF90(grid%xmol, x_p, ierr)
  call VecGetArrayF90(grid%sat, s_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecGetArrayF90(grid%h, h_p, ierr)
  call VecGetArrayF90(grid%u, u_p, ierr)
  call VecGetArrayF90(grid%pcw, pc_p, ierr)
  call VecGetArrayF90(grid%kvr, kvr_p, ierr)
  call VecGetArrayF90(grid%hen, hen_p,ierr) 
  call VecGetArrayF90(grid%df, df_p,ierr) 
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)
! print *,'initadj gotten pointers' 


 do n = 1, grid%nlmax
        jn = 1 + (n-1)*grid%nphase
        ii1=1+(n-1)*grid%nphase; ii2=n*grid%nphase
        iicap=int(icap_p(n))
        i1_hencoeff=2+(n-1)*grid%nphase*grid%nspec
        p_p(jn)=p_p(jn+1)
        x_p(jn)=x_p(jn+1)        ! For debug
        s_p(jn)=1.D0-s_p(jn+1)   

        iphase_p(n)=6
        if(s_p(jn+1)<eps) iphase_p(n)=2
        if(s_p(jn+1)>1.D0-eps) iphase_p(n)=4
         
        density_p(ii2)=p_p(jn+1) / (t_p(n)+273.15D0) / Rg / 1D3
      
!       print *,'init:: ',n,jn, p_p(jn+1),t_p(n),s_p(jn+1), density_p(ii2)
       
       
        call mixture_eos_noderiv (p_p(jn+1),t_p(n),1.D0, s_p(jn+1),&
          grid%scale,grid%nphase,grid%nspec,grid%npricomp,&
          iicap, grid%swir(iicap), grid%lambda(iicap),&
          grid%alpha(iicap),grid%pckrm(iicap),&
          grid%pcwmax(iicap),grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
          density_p(ii1:ii2),avgmw_p(ii1:ii2),& 
          h_p(ii1:ii2),u_p(ii1:ii2),&
          df_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
          hen_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
          grid%xxphi_co2(n),&
          pc_p(ii1:ii2),kvr_p(ii1:ii2),&
          ierr,grid%itable)

!    print *,'initadj',eps,sat_pressure,jn  
! 1: Assume initially the reservoir is under thermal equilibrium
     
      if((s_p(jn+1)>eps).and.(s_p(jn+1)<(1.D0-eps)))then  
!       if((s_p(jn+1)<(1.D0-eps)))then  
        ! xga=(p_p(jn)-sat_pressure)/(p_p(jn)-sat_pressure*hen_p(i1_hencoeff))
         xga=1.D0-grid%yh2o_in_co2 ! for super-critical co2 phase
         xgw=1.D0-xga
   
         xgw=1.D0-xga
        ! if(xga<x_p(jn+1)) xga=x_p(jn+1)

         xla=xga * hen_p(i1_hencoeff)
         xlw=1.D0-xla
         
 
         if(abs(xgw+xga-1.D0)>eps) then
            print *,'Wrong assignment in init condition, STOP'
            stop
         end if
         
         x_p(jn) = xla
!         x_p(i1_hencoeff    ) = xla
         x_p(jn+1 ) = xga
!         x_p(i1_hencoeff + 2) = xga
! remember xmol is nphase degree of freedom, need to be changed later 
!     print *,'initadj_2ph',jn,eps,sat_pressure,xga 
      call mixture_eos_noderiv (p_p(jn+1),t_p(n), xga , s_p(jn+1),&
        grid%scale,grid%nphase,grid%nspec,grid%npricomp,&
        iicap, grid%swir(iicap), grid%lambda(iicap),&
        grid%alpha(iicap),grid%pckrm(iicap),&
        grid%pcwmax(iicap),grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
        density_p(ii1:ii2),avgmw_p(ii1:ii2),& 
        h_p(ii1:ii2),u_p(ii1:ii2),&
        df_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        hen_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
        grid%xxphi_co2(n),&
        pc_p(ii1:ii2),kvr_p(ii1:ii2),&
        ierr,grid%itable)

!       print *,'initend',p_p(jn),sat_pressure,xgw,xlw,xga,xla,hen_p(i1_hencoeff)

      end if
      

     ! if(n==1)then
         !print *, ' init adjustmet::__________________________'
         !print *,n,ii1,jn,iicap, p_p(jn+1),t_p(n), s_p(jn+1)
         !print *, x_p(i1_hencoeff - 1:i1_hencoeff + 2)
         !print *,  grid%scale,grid%nphase,grid%nspec,grid%npricomp
         !print *,'Cap  ',  iicap, grid%swir(iicap), grid%lambda(iicap),&
!                       grid%alpha(iicap),grid%pckrm(iicap),&
!                       grid%pcwmax(iicap)
        !print *, 'Ps  ', sat_pressure
        !print *, 'den   ', density_p(ii1:ii2),avgmw_p(ii1:ii2)
        !print *, 'H  U    ',   h_p(ii1:ii2),u_p(ii1:ii2)
        !print *,'diff   ', df_p(1+(n-1)*grid%nphase*grid%nspec:&
!                            n*grid%nphase*grid%nspec)
        !print *,'Henry  ',hen_p(1+(n-1)*grid%nphase*grid%nspec:&
!                            n*grid%nphase*grid%nspec)
        !print *,'pc  ', pc_p(ii1:ii2)
        !print *,'kvr   ', kvr_p(ii1:ii2)
    !  endif

 enddo

   do nc = 1, grid%nconnbc

       m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
       

       if(m<0)then
          print *, "Wrong boundary node index... STOP!!!"
          stop
       end if

       ibc = grid%ibconn(nc)
       
!      print *,'initadj_bc',nc,ibc,grid%ibndtyp(ibc),grid%nconnbc

       if(grid%ibndtyp(ibc)==1)then
          iicap=int(icap_p(m))


          grid%density_bc(2)=grid%pressurebc(2,nc) /(grid%tempbc(nc)+273.15D0) / Rg * 1.D-3
          
!         print *,'pflow_2ph_mix_noder: ',ibc,grid%density_bc,rg
          
          call mixture_eos_noderiv (grid%pressurebc(2,nc),grid%tempbc(nc),&
            1.0D0,grid%sgbc(nc),&
            grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
            iicap, grid%swir(iicap),grid%lambda(iicap),&
            grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), &
            grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
            grid%density_bc,grid%avgmw_bc,grid%hh_bc,grid%uu_bc,&
            grid%df_bc,grid%hen_bc,grid%xxphi_co2_bc(nc),grid%pc_bc,grid%kvr_bc,&
            ierr,grid%itable)
  
!         print *,'initadj_bc',nc,ibc,grid%ibndtyp(ibc), sat_pressure,jn
          
          if((grid%sgbc(nc)>eps).and.(grid%sgbc(nc)<(1.D0-eps)))then  
 !          if((grid%sgbc(ibc)<(1.D0-eps)))then  
            xga=(grid%pressurebc(2,nc)-sat_pressure)/ &
                 (grid%pressurebc(2,nc)-sat_pressure*grid%hen_bc(2))
            
             xgw=1.D0-xga
             ! if(xga<x_p(jn+1)) xga=x_p(jn+1)

             xla=xga * grid%hen_bc(2)
             xlw=1.D0-xla
             
 
             if(abs(xgw+xga-1.D0)>eps) then
                print *,'Wrong assignment in init condition, STOP'
                stop
             end if

             
             grid%concbc(nc) = xga

!           print *,'pflow_2pha_noderiv: ',ibc,grid%density_bc
          
            call mixture_eos_noderiv (grid%pressurebc(2,nc),grid%tempbc(nc),&
              grid%concbc(nc),grid%sgbc(nc),&
              grid%scale,grid%nphase,grid%nspec,grid%npricomp, &
              iicap, grid%swir(iicap),grid%lambda(iicap),&
              grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap), &
              grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
              grid%density_bc,grid%avgmw_bc,grid%hh_bc,grid%uu_bc,&
              grid%df_bc,grid%hen_bc,grid%xxphi_co2_bc(nc),grid%pc_bc,grid%kvr_bc,&
              ierr,grid%itable)
          end if
      end if

    enddo

  call VecRestoreArrayF90(grid%pressure, p_p, ierr)
  call VecRestoreArrayF90(grid%temp, t_p, ierr)
  call VecRestoreArrayF90(grid%xmol, x_p, ierr)
  call VecRestoreArrayF90(grid%sat, s_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecRestoreArrayF90(grid%h, h_p, ierr)
  call VecRestoreArrayF90(grid%u, u_p, ierr)
  call VecRestoreArrayF90(grid%pcw, pc_p, ierr)
  call VecRestoreArrayF90(grid%kvr, kvr_p, ierr)
  call VecRestoreArrayF90(grid%hen, hen_p,ierr) 
  call VecRestoreArrayF90(grid%df, df_p,ierr) 
  call VecRestoreArrayF90(grid%icap, icap_p, ierr)
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  
  call VecCopy(grid%density,grid%ddensity,ierr)
  call VecCopy(grid%iphas,grid%iphas_old,ierr)
   
end subroutine pflow_2phase_initadj


 subroutine pflow_2phase_massbal(grid)
 
  use mixture_module  
  implicit none
  type(pflowGrid) :: grid 
  
 
  integer :: ierr,icall
  integer :: n
  integer :: j, jn
  integer :: ii1,ii2,iicap
 
#include "definitions.h"

  PetscScalar, pointer :: yy_p(:)
  
  PetscScalar, pointer ::  porosity_p(:), volume_p(:), &
                           density_p(:), avgmw_p(:), u_p(:), &
                           h_p(:), df_p(:), hen_p(:),&
                           pc_p(:), kvr_p(:),&
                           ithrm_p(:),icap_p(:),iphase_p(:) 

!  integer, pointer ::iphase_p(:)
  
  real*8 :: sat_pressure, pvol, satw  ! Saturation pressure of water.
 
  real*8 :: totl,totg,totl0,totg0, accl,accg, xxlw,xxgw,xxla,xxga
  
  data icall/0/

  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%porosity, porosity_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecGetArrayF90(grid%h, h_p, ierr)
  call VecGetArrayF90(grid%u, u_p, ierr)
  call VecGetArrayF90(grid%pcw, pc_p, ierr)
  call VecGetArrayF90(grid%kvr, kvr_p, ierr)
  call VecGetArrayF90(grid%hen, hen_p,ierr) 
  call VecGetArrayF90(grid%df, df_p,ierr) 
  call VecGetArrayF90(grid%icap, icap_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%iphas, iphase_p, ierr)

  totl=0.d0
  totg=0.D0
  do n = 1, grid%nlmax
    jn = 1 + (n-1)*grid%nphase
    ii1=1+(n-1)*grid%nphase; ii2=n*grid%nphase
    iicap=int(icap_p(n))

    call mixture_eos_noderiv (PRESSURE(n),TEMP(n),CONC(n),SATG(n),&
      grid%scale,grid%nphase,grid%nspec,grid%npricomp,&
      iicap, grid%swir(iicap),grid%lambda(iicap),&
      grid%alpha(iicap),grid%pckrm(iicap),&
      grid%pcwmax(iicap),grid%pcbetac(iicap),grid%pwrprm(iicap),sat_pressure,&
      density_p(ii1:ii2),avgmw_p(ii1:ii2),& 
      h_p(ii1:ii2),u_p(ii1:ii2),&
      df_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
      hen_p(1+(n-1)*grid%nphase*grid%nspec:n*grid%nphase*grid%nspec),&
      grid%xxphi_co2(n),&
      pc_p(ii1:ii2),kvr_p(ii1:ii2),&
      ierr,grid%itable)
  enddo


!---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
  ! ng = grid%nL2G(n)   ! corresponding ghost index
   
    pvol=volume_p(n)*porosity_p(n)
    SATW=1.D0-SATG(n)

!   Pressure equation accumulation term
    accl = 0.d0
    accg = 0.d0
!   do j = 1, grid%nphase
    j=1
    jn = j + (n-1)*grid%nphase
  
    xxga = CONC(n)
    xxgw = 1.D0-xxga
    xxla =  xxga*hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec)
    xxlw = 1.D0 - xxla
!    yylw = CCONC(n)
!    yyla = 1.D0-yylw
!    yyga = henry_p(n) * yyla
!    yygw = 1.D0 - yyga
!    print *,'initaccum', xxla,hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec)
    accl=pvol * density_p(jn) * xxla * SATW
    accg=pvol * density_p(jn+1) * xxga * SATG(n)
    totl=totl+accl
    totg=totg+accg
    
!   print *,'totco2: ',n,jn,xxga,xxla,pvol,satw,density_p(jn),density_p(jn+1)
    
 !print *,n,accum_p(p1),accum_p(t1),accum_p(c1),accum_p(s1)
 !print *,  n, PRESSURE(n),TEMP(n), density_p(jn), density_p(jn+1), u_p(jn),u_p(jn+1),&
 !hen_p(2+(j-1)*grid%nspec+(n-1)*grid%nphase*grid%nspec),kvr_p(jn),kvr_p(jn+1)

  end do

  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%porosity, porosity_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%avgmw, avgmw_p, ierr)
  call VecRestoreArrayF90(grid%h, h_p, ierr)
  call VecRestoreArrayF90(grid%u, u_p, ierr)
  call VecRestoreArrayF90(grid%pcw, pc_p, ierr)
  call VecRestoreArrayF90(grid%kvr, kvr_p, ierr)
  call VecRestoreArrayF90(grid%hen, hen_p,ierr) 
  call VecRestoreArrayF90(grid%df, df_p,ierr)
  call VecRestoreArrayF90(grid%icap, icap_p, ierr) 
  call VecRestoreArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p, ierr)
  
  if(grid%commsize >1)then
    call MPI_REDUCE(totl, totl0,1, MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_BCAST(totl0,1, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(totg, totg0,1, MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_BCAST(totg0,1, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    totl = totl0
    totg = totg0
  endif 
  
  
  if(grid%myrank==0)then
    write(*,'(" Total CO2: t=",1pe13.6,'' liq:'',1pe13.6,'' gas:'',1pe13.6, &
 &  '' tot:'',1pe13.6,'' [kmol]'')') grid%t/grid%tconv,totl,totg,totl+totg
    write(IUNIT2,'(" Total CO2: t=",1pe13.6,'' liq:'',1pe13.6,'' gas:'',1pe13.6, &
 &  '' tot:'',1pe13.6,'' [kmol]'')') grid%t/grid%tconv,totl,totg,totl+totg
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,*) '# time   dt   totl   totg   tot'
      icall = 1
    endif
    write(13,'(1p5e12.4)') grid%t/grid%tconv,grid%dt/grid%tconv,totl,totg,totl+totg
  endif   
end subroutine pflow_2phase_massbal




end module TTPHASE_module

#undef PPRESSURE_LOC
#undef PPRESSURE      
#undef PRESSURE       
#undef TTEMP_LOC
#undef TTEMP
#undef TEMP
#undef CCONC_LOC
#undef CCONC
#undef CONC
#undef SSATG_LOC
#undef SSATG
#undef SATG
