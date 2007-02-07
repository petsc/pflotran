!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_init.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_init.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:33:01  lichtner
! Revised storage of phik and surf. Set rho=1 temporarily.
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module ptran_init_module

  public

  interface ptran_init
    module procedure ptran_init
  end interface ptran_init

contains

  subroutine ptran_init (da,da_mat,da_1dof,da_kin,sles,ksp)

  use ptran_global_module
  use trdynmem_module
  use water_eos_module

  implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#ifdef USE_PETSC216
    ! petsc-2.1.6
#include "include/finclude/petscsles.h"
#endif
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA :: da, da_mat, da_1dof, da_kin

#ifdef USE_PETSC216
  SLES :: sles
#endif

#ifdef USE_PETSC221
  integer :: sles
#endif

! PetscViewer   vecviewer

  KSPType :: ksp_type
  PCType  :: pc_type
  KSP   ::  ksp
  PC    ::  pc
  Vec :: temp1_nat_vec,temp2_nat_vec,temp3_nat_vec,temp4_nat_vec
  Vec :: phik_tmp, surf_tmp

  integer :: ierr,i,j,k,ir,n,ng,nr, na
! integer :: isumjl
! integer :: dfill(ncomp*ncomp), ofill(ncomp*ncomp)
  
  real*8, pointer :: temploc_p(:),pressloc_p(:)
  real*8 :: val
  
  real*8, pointer :: dx_p(:), dy_p(:), dz_p(:)
  
  ierr = 0

!--convert units
! keff = keff/ceq*1.d3 ! convert to mol/L/s
  difaq = difaq*1.d-4  ! cm^2/s -> m^2/s
  difgas = difgas*1.d-4! cm^2/s -> m^2/s
  vlx0 = vlx0/yrsec    ! m/y -> m/s
  vly0 = vly0/yrsec    ! m/y -> m/s
  vlz0 = vlz0/yrsec    ! m/y -> m/s
  vgx0 = vgx0/yrsec    ! m/y -> m/s
  vgy0 = vgy0/yrsec    ! m/y -> m/s
  vgz0 = vgz0/yrsec    ! m/y -> m/s

! do i = 1, kplot
!   tplot(i) = tplot(i)*yrsec ! y -> s
! enddo
! dt    = dt*yrsec     ! y -> s
! dtmax = dtmax*yrsec  ! y -> s

  nmat = ncomp

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,1,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       da_1dof,ierr)
  
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,nmat,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       da_mat,ierr)
  
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,ncomp,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       da,ierr)

  if (nkin > 0) &
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,nkin,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       da_kin,ierr)

!---petsc-2.1.6
  if (iblkfmt == 1) then ! blocked
!   call DAGetMatrix(da,MATBAIJ,A,ierr)
    call DAGetMatrix(da_mat,MATBAIJ,A,ierr)
  else if (iblkfmt == 0) then ! compressed storage
!   call DAGetMatrix(da,MATAIJ,A,ierr)
    call DAGetMatrix(da_mat,MATAIJ,A,ierr)
  else if (iblkfmt == 2) then ! element storage
    call DAGetMatrix(da,MATAIJ,A,ierr)
  endif

!===============================================================
! natural vector:  physical space (natural) ordering
! global vector:  local processor ordering
! local vector:   local/ghosted processor ordering (stored as sequential vec)

! global -> local
! call DAGlobalToLocalBegin(da,cc,INSERT_VALUES,ccloc,ierr)
! call DAGlobalToLocalEnd(da,cc,INSERT_VALUES,ccloc,ierr)

! VecGetArrayF90: global -> local without ghosted nodes
! VecGetArrayF90: local -> local with ghosted nodes
!===============================================================

  call DACreateGlobalVector(da_mat,b,ierr)
! call DACreateGlobalVector(da,b,ierr)
  call VecDuplicate(b,x,ierr)

  call DACreateGlobalVector(da,c,ierr)
  call VecDuplicate(c,cc,ierr)
  call VecDuplicate(c,dc,ierr)
  call VecDuplicate(c,ptran_c0,ierr)
  
  call DACreateGlobalVector(da_1dof,porosity,ierr)
  call VecDuplicate(porosity,por,ierr)
  call VecDuplicate(porosity,tortuosity,ierr)
  call VecDuplicate(porosity,temp,ierr)
  call VecDuplicate(porosity,press,ierr)
  call VecDuplicate(porosity,sat,ierr)
  call VecDuplicate(porosity,ssat,ierr)
  call VecDuplicate(porosity,dx,ierr)
  call VecDuplicate(porosity,dy,ierr)
  call VecDuplicate(porosity,dz,ierr)
  call VecDuplicate(porosity,nreg,ierr)
  call VecDuplicate(porosity,ghost,ierr)
  call VecDuplicate(porosity,xphi_co2,ierr)
  call VecDuplicate(porosity,den_co2,ierr)
  
  if (nkin > 0) then
    call DACreateGlobalVector(da_kin,phik,ierr)
    call VecDuplicate(phik,surf,ierr)
    call VecDuplicate(phik,phik0,ierr)
    call VecDuplicate(phik,surf0,ierr)
    call VecDuplicate(phik,rkin,ierr)
    call VecDuplicate(phik,rrkin,ierr)
  endif

  call DACreateLocalVector(da,ccloc,ierr)
  
  call DACreateLocalVector(da_1dof,porloc,ierr)
  call VecDuplicate(porloc,tort_loc,ierr)
  call VecDuplicate(porloc,temploc,ierr)
  call VecDuplicate(porloc,pressloc,ierr)
  call VecDuplicate(porloc,sat_loc,ierr)
  call VecDuplicate(porloc,ssat_loc,ierr)
  call VecDuplicate(porloc,xphi_co2_loc,ierr)
  call VecDuplicate(porloc,den_co2_loc,ierr)
  call VecDuplicate(porloc,dx_loc,ierr)
  call VecDuplicate(porloc,dy_loc,ierr)
  call VecDuplicate(porloc,dz_loc,ierr)
  call VecDuplicate(porloc,ghost_loc,ierr)

#ifdef USE_PETSC216
    ! petsc-2.1.6
  call SLESCreate(PETSC_COMM_WORLD,sles,ierr)
#endif

#ifdef USE_PETSC221
    ! petsc-2.2.0
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
#endif

! call DAGetCorners(da,nxs,nys,nzs,nlx,nly,nlz,ierr)
  call DAGetCorners(da_mat,nxs,nys,nzs,nlx,nly,nlz,ierr)
!-----------------------------------------------------------------
! nxs,nys,nzs = global corner (starting) indices (without ghost nodes)
! nlx,nly,nlz = local grid dimension (without ghost nodes)
!-----------------------------------------------------------------

  nxe   = nxs+nlx
  nye   = nys+nly
  nze   = nzs+nlz
  nlxy  = nlx*nly
  nlxz  = nlx*nlz
  nlyz  = nly*nlz
  nlmax = nlx*nly*nlz
  nldof = nlmax*ncomp

! Note that we need to allocate this array, even when running coupled 
! with pflow, because the flow code uses a PETSc vector to store the 
! volumes.
! if (using_pflowGrid == PETSC_FALSE) &
  allocate(vb(nlmax))

  call DAGetGhostCorners(da_mat,ngxs,ngys,ngzs,ngx,ngy,ngz,ierr)
! call DAGetGhostCorners(da,ngxs,ngys,ngzs,ngx,ngy,ngz,ierr)
!-----------------------------------------------------------------
! ngxs,ngys,ngzs = global corner (starting) indices (with ghost nodes)
! ngx,ngy,ngz = local grid dimension (with ghost nodes)
!-----------------------------------------------------------------

  ngxe  = ngxs+ngx
  ngye  = ngys+ngy
  ngze  = ngzs+ngz
  ngxy  = ngx*ngy
  ngxz  = ngx*ngz
  ngyz  = ngy*ngz
  ngmax = ngx*ngy*ngz
  ngdof = ngmax*ncomp
  
  write(*,'(" myrank= ",i3,", nlmax= ",i5,", nlx,y,z= ",3i4, &
  & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
  myrank,nlmax,nlx,nly,nlz,nxs,nxe,nys,nye,nzs,nze

  write(*,'(" myrank= ",i3,", ngmax= ",i5,", ngx,y,z= ",3i4, &
  & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
  myrank,ngmax,ngx,ngy,ngz,ngxs,ngxe,ngys,ngye,ngzs,ngze

! Compute arrays for indexing between local ghosted and non-ghosted 
! arrays
! if (using_pflowGrid .ne. PETSC_TRUE) then
    allocate(nL2A(nlmax))
	allocate(nL2G(nlmax))
    allocate(nG2L(ngmax))

    nG2L(1:ngmax) = 0

    istart = nxs-ngxs
    jstart = nys-ngys
    kstart = nzs-ngzs
    iend = istart+nlx-1
    jend = jstart+nly-1
    kend = kstart+nlz-1

  ! Local <-> Ghosted Transformation
    n = 0
    do k=kstart,kend
      do j=jstart,jend
        do i=istart,iend
          n = n + 1
          ng = i+j*ngx+k*ngxy+1
!         vb(n) = dx0(i)*dy0(j)*dz0(k)
          nL2G(n) = ng
          nG2L(ng) = n
        enddo
      enddo
    enddo

    do i=1,ngmax
      j = nG2L(i)
      if (j > 0) then
        k = nL2G(j)
        if(i /= k) then
          print *,'Error in ghost-local node numbering for ghost node =', i
          print *,'node_id_gtol(i) =', j
          print *,'node_id_ltog(node_id_gtol(i)) =', k
          stop
        endif
      endif
    enddo
	
	 !Local(non ghosted)->Natural(natural order starts from 0)
    n=0
    do k=1,nlz
      do j=1,nly
        do i=1,nlx
          n = n + 1
          na = i-1+nxs+(j-1+nys)*nx+(k-1+nzs)*nxy
          if(na>(nmax-1)) print *,'Wrong Nature order....'
          nL2A(n) = na
          !print *,grid%myrank, k,j,i,n,na
          !grid%nG2N(ng) = na
        enddo
      enddo
    enddo

	
	
	
! endif
  
  ! identify corner ghost nodes
  call VecSet(ghost_loc,-1.d0,ierr)
  call VecSet(ghost,0.d0,ierr)
  call DAGlobalToLocalBegin(da_1dof,ghost,INSERT_VALUES,ghost_loc,ierr)
  call DAGlobalToLocalEnd(da_1dof,ghost,INSERT_VALUES,ghost_loc,ierr)
  call VecDestroy(ghost,ierr)
  call VecGetArrayF90(ghost_loc,ghost_loc_p,ierr)
! do n = 1, ngmax
!   m = nG2L(n)
!   print *,'ptraninit-ghost: ',myrank,n,m,ghost_loc_p(n)
! enddo
! call VecRestoreArrayF90(ghost_loc,ghost_loc_p,ierr)
  
!-------preconditioner-------------------------
#ifdef USE_PETSC216
  ! petsc-2.1.6
  call SLESGetPC(sles, pc, ierr)
#endif

#ifdef USE_PETSC221
  ! petsc-2.2.0
  call KSPGetPC(ksp, pc, ierr)
#endif

! pc_type = PCILU
  if (iblkfmt == 0) then
    pc_type = PCJACOBI
  else
    pc_type = PCBJACOBI
  endif
! pc_type = PCSLES
! pc_type = PCASM
! pc_type = PCNONE
  call PCSetType(pc,pc_type,ierr)

! call PetscOptionsSetValue('-pc_ilu_damping','1.d-10',ierr)
! call PCILUSetDamping(pc,1.d-14,ierr)

!-------krylov subspace method ----------------
#ifdef USE_PETSC221
  call KSPSetFromOptions(ksp,ierr)
#endif

#ifdef USE_PETSC216
    ! petsc-2.1.6
  call SLESSetFromOptions(sles,ierr)
  call SLESGetKSP(sles,ksp,ierr)
#endif

  call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

! ksp_type = KSPGMRES
  ksp_type = KSPFGMRES
  call KSPSetType(ksp,ksp_type,ierr)

! Solver
! atol_petsc = 1.d-12
! dtol_petsc = 1.d+5
! rtol_petsc = 1.d-12
! maxits_petsc = 100
! restart_petsc = 30

! KSPDefaultConverged() reaches convergence when
!    rnorm < MAX (rtol * rnorm_0, atol);
! Divergence is detected if
!    rnorm > dtol * rnorm_0,


  call KSPSetTolerances(ksp,rtol_petsc,atol_petsc,dtol_petsc, &
      maxits_petsc,ierr)

!--set initial conditions

  if (myrank == 0) write(*,*) '--> allocate memory'
  call trallocate
  
! set dx, dy, dz, nreg
  call DACreateNaturalVector(da_1dof,temp1_nat_vec,ierr)
  call VecDuplicate(temp1_nat_vec, temp2_nat_vec, ierr)
  call VecDuplicate(temp1_nat_vec, temp3_nat_vec, ierr)
  call VecDuplicate(temp1_nat_vec, temp4_nat_vec, ierr)
  if (myrank == 0) then
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          n = i+(j-1)*nx+(k-1)*nxy-1
          val = dx0(i)
          call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
          val = dy0(j)
          call VecSetValue(temp2_nat_vec,n,val,INSERT_VALUES,ierr)
          val = dz0(k)
          call VecSetValue(temp3_nat_vec,n,val,INSERT_VALUES,ierr)
          val = nreg_val(n+1)
          call VecSetValue(temp4_nat_vec,n,val,INSERT_VALUES,ierr)
        enddo
      enddo
    enddo
  endif

  call VecAssemblyBegin(temp1_nat_vec,ierr)
  call VecAssemblyEnd(temp1_nat_vec,ierr)
  call VecAssemblyBegin(temp2_nat_vec,ierr)
  call VecAssemblyEnd(temp2_nat_vec,ierr)
  call VecAssemblyBegin(temp3_nat_vec,ierr)
  call VecAssemblyEnd(temp3_nat_vec,ierr)
  call VecAssemblyBegin(temp4_nat_vec,ierr)
  call VecAssemblyEnd(temp4_nat_vec,ierr)
  
  call DANaturalToGlobalBegin(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                                dx,ierr)
  call DANaturalToGlobalEnd(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                              dx,ierr)
  call DANaturalToGlobalBegin(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                                dy,ierr)
  call DANaturalToGlobalEnd(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                              dy,ierr)
  call DANaturalToGlobalBegin(da_1dof,temp3_nat_vec,INSERT_VALUES, &
                                dz,ierr)
  call DANaturalToGlobalEnd(da_1dof,temp3_nat_vec,INSERT_VALUES, &
                              dz,ierr)
  call DANaturalToGlobalBegin(da_1dof,temp4_nat_vec,INSERT_VALUES, &
                                nreg,ierr)
  call DANaturalToGlobalEnd(da_1dof,temp4_nat_vec,INSERT_VALUES, &
                              nreg,ierr)
  
  call VecDestroy(temp1_nat_vec,ierr)
  call VecDestroy(temp2_nat_vec,ierr)
  call VecDestroy(temp3_nat_vec,ierr)
  call VecDestroy(temp4_nat_vec,ierr)
      
  
  ! Calculate volumes of local cells.
  if (igeom == 1) then
    call VecGetArrayF90(dx, dx_p, ierr)
    call VecGetArrayF90(dy, dy_p, ierr)
    call VecGetArrayF90(dz, dz_p, ierr)
    do n = 1, nlmax
      vb(n) = dx_p(n) * dy_p(n) * dz_p(n)
    enddo
    call VecRestoreArrayF90(dx, dx_p, ierr)
    call VecRestoreArrayF90(dy, dy_p, ierr)
    call VecRestoreArrayF90(dz, dz_p, ierr)
  else if (igeom == 2) then
    call VecGetArrayF90(dz, dz_p, ierr)
    allocate(rd(0:nx))
    rd = 0.d0
    rd(0) = 0.d0 
    do i = 1, nx
      rd(i) = rd(i-1) + dx0(i)
    enddo
    do n=1, nlmax
      i = mod(mod((n),nlxy),nlx)!+(grid%ngxs-grid%nxs)
      if (i==0) i = nlx
      vb(n) = Pi * (rd(i+nxs) + rd(i-1+nxs))*&
      (rd(i+nxs) - rd(i-1+nxs)) * dz_p(n)
 !     print *, 'setup: Vol ', myrank, n,i, rd(i+nxs),vb(n)
    enddo
    call VecRestoreArrayF90(dz, dz_p, ierr)
  else if (igeom == 3) then
    allocate(rd(0:nx))
    rd = 0.d0
    rd(0) = 0.d0 
    do i = 1, nx
      rd(i) = rd(i-1) + dx0(i)
    enddo
    do i = 1, nlmax
      vb(n) = 4.d0*Pi*(rd(i)+rd(i-1))*(rd(i)-rd(i-1))/3.d0
    enddo
  endif
  
! set temperature, pressure, saturation, and porosity    
  if (iregfld == 0) then

    if (using_pflowGrid .ne. PETSC_TRUE) then
      call VecSet(porosity,por0,ierr)
      call VecSet(temp,temp0+tkelvin,ierr)
    endif
    
    call VecSet(press,pref0,ierr)
    call VecSet(sat,sat0,ierr)
    call VecSet(ssat,sat0,ierr)
    
    call VecSet(tortuosity,tor0,ierr)

    call density(temp0,pref0,rho0)
    rho0 = rho0*1.d-3
    rho0 = 1.d0

    do n = 1, ngmax
      if (ghost_loc_p(n) == -1) cycle
      rho(n) = rho0
    enddo
    
  else if (iregfld > 0) then
  
    call VecSet(sat,sat0,ierr)
    call VecSet(ssat,sat0,ierr)

!---set porosity, pressure, temperature, and saturation by region
    call DACreateNaturalVector(da_1dof,temp1_nat_vec,ierr)
    call VecDuplicate(temp1_nat_vec, temp2_nat_vec, ierr)
    call VecDuplicate(temp1_nat_vec, temp3_nat_vec, ierr)
    call VecDuplicate(temp1_nat_vec, temp4_nat_vec, ierr)
    if (myrank == 0) then
      do ir = 1,iregfld
        do k = k1reg(ir),k2reg(ir)
          do j = j1reg(ir),j2reg(ir)
            do i = i1reg(ir),i2reg(ir)
              n = i+(j-1)*nx+(k-1)*nxy-1
              val = por_reg(ir)
              call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
              val = tor_reg(ir)
              call VecSetValue(temp2_nat_vec,n,val,INSERT_VALUES,ierr)
              val = pref_reg(ir)
              call VecSetValue(temp3_nat_vec,n,val,INSERT_VALUES,ierr)
              val = temp_reg(ir)+tkelvin
              call VecSetValue(temp4_nat_vec,n,val,INSERT_VALUES,ierr)
            enddo
          enddo
        enddo
      enddo
    endif
    
    call VecAssemblyBegin(temp1_nat_vec,ierr)
    call VecAssemblyEnd(temp1_nat_vec,ierr)

    call VecAssemblyBegin(temp2_nat_vec,ierr)
    call VecAssemblyEnd(temp2_nat_vec,ierr)
  
    call VecAssemblyBegin(temp3_nat_vec,ierr)
    call VecAssemblyEnd(temp3_nat_vec,ierr)
  
    call VecAssemblyBegin(temp4_nat_vec,ierr)
    call VecAssemblyEnd(temp4_nat_vec,ierr)
  
    call DANaturalToGlobalBegin(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                                porosity,ierr)
    call DANaturalToGlobalEnd(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                              porosity,ierr)
    call DANaturalToGlobalBegin(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                                tortuosity,ierr)
    call DANaturalToGlobalEnd(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                              tortuosity,ierr)
    call DANaturalToGlobalBegin(da_1dof,temp3_nat_vec,INSERT_VALUES, &
                                press,ierr)
    call DANaturalToGlobalEnd(da_1dof,temp3_nat_vec,INSERT_VALUES, &
                              press,ierr)
    call DANaturalToGlobalBegin(da_1dof,temp4_nat_vec,INSERT_VALUES, &
                                temp,ierr)
    call DANaturalToGlobalEnd(da_1dof,temp4_nat_vec,INSERT_VALUES, &
                              temp,ierr)
    
    call VecDestroy(temp1_nat_vec,ierr)
    call VecDestroy(temp2_nat_vec,ierr)
    call VecDestroy(temp3_nat_vec,ierr)
    call VecDestroy(temp4_nat_vec,ierr)
    
!   call VecView(temp,PETSC_VIEWER_STDOUT_WORLD,ierr)

!---compute density
    call DAGlobalToLocalBegin(da_1dof,temp,INSERT_VALUES,temploc,ierr)
    call DAGlobalToLocalEnd(da_1dof,temp,INSERT_VALUES,temploc,ierr)
    
    call DAGlobalToLocalBegin(da_1dof,press,INSERT_VALUES,pressloc,ierr)
    call DAGlobalToLocalEnd(da_1dof,press,INSERT_VALUES,pressloc,ierr)
  
    call VecGetArrayF90(temploc,temploc_p,ierr)
    call VecGetArrayF90(pressloc,pressloc_p,ierr)
    do n = 1, ngmax
      if (ghost_loc_p(n) == -1) cycle
      call density (temploc_p(n),pressloc_p(n),rho0)
!     rho(n) = rho0*1.d-3
      rho(n) = 1.d0
    enddo
    call VecRestoreArrayF90(temploc,temploc_p,ierr)
    call VecRestoreArrayF90(pressloc,pressloc_p,ierr)
  endif
  
! set mineral volume fraction and area by region
  if (nkin > 0) then
!   rkin = zero
    call VecSet(rkin,zero,ierr)
    call DACreateGlobalVector(da_1dof, phik_tmp, ierr)
    call VecDuplicate(phik_tmp, surf_tmp, ierr)
    call DACreateNaturalVector(da_1dof,temp1_nat_vec,ierr)
    call VecDuplicate(temp1_nat_vec,temp2_nat_vec,ierr)
    do nr = 1, nkin
      if (myrank == 0) then
        do ir = 1,iregkin(nr)
          do k = k1kin(nr,ir),k2kin(nr,ir)
            do j = j1kin(nr,ir),j2kin(nr,ir)
              do i = i1kin(nr,ir),i2kin(nr,ir)
                n = i+(j-1)*nx+(k-1)*nxy - 1
                val = phik_reg(nr,ir)
                call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
                val = surf_reg(nr,ir)
                call VecSetValue(temp2_nat_vec,n,val,INSERT_VALUES,ierr)
              enddo
            enddo
          enddo
        enddo
      endif

      call VecAssemblyBegin(temp1_nat_vec,ierr)
      call VecAssemblyEnd(temp1_nat_vec,ierr)

      call VecAssemblyBegin(temp2_nat_vec,ierr)
      call VecAssemblyEnd(temp2_nat_vec,ierr)

      call DANaturalToGlobalBegin(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                                phik_tmp,ierr)
      call DANaturalToGlobalEnd(da_1dof,temp1_nat_vec,INSERT_VALUES, &
                              phik_tmp,ierr)

      call DANaturalToGlobalBegin(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                                surf_tmp,ierr)
      call DANaturalToGlobalEnd(da_1dof,temp2_nat_vec,INSERT_VALUES, &
                              surf_tmp,ierr)

      call VecStrideScatter(phik_tmp,nr-1,phik,INSERT_VALUES,ierr)
      call VecStrideScatter(surf_tmp,nr-1,surf,INSERT_VALUES,ierr)
    enddo
    
    call VecDestroy(temp1_nat_vec,ierr)
    call VecDestroy(temp2_nat_vec,ierr)
    call VecDestroy(phik_tmp,ierr)
    call VecDestroy(surf_tmp,ierr)
    
    call VecCopy(phik,phik0,ierr)
    call VecCopy(surf,surf0,ierr)

!   deallocate(surf_reg)

#if 0
    do nr = 1, nkin
      do ir = 1,iregkin(nr)
      
        kk1 = k1kin(nr,ir) - nzs
        kk2 = k2kin(nr,ir) - nzs
        jj1 = j1kin(nr,ir) - nys
        jj2 = j2kin(nr,ir) - nys
        ii1 = i1kin(nr,ir) - nxs
        ii2 = i2kin(nr,ir) - nxs
        
        kk1 = max(1,kk1)
        kk2 = min(nlz,kk2)
        jj1 = max(1,jj1)
        jj2 = min(nly,jj2)
        ii1 = max(1,ii1)
        ii2 = min(nlx,ii2)
        
        if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
        
        do k = kk1,kk2
          do j = jj1,jj2
            do i = ii1,ii2
              n = i+(j-1)*nlx+(k-1)*nlxy
!             phik(nr,n) = phik_reg(nr,ir)
              surf(nr,n) = surf_reg(nr,ir)
            enddo
          enddo
        enddo
      enddo
    enddo
#endif

    call VecGetArrayF90(phik,phik_p,ierr)
    call VecGetArrayF90(surf,surf_p,ierr)
    call VecGetArrayF90(phik0,phik0_p,ierr)
    call VecGetArrayF90(surf0,surf0_p,ierr)
    call VecGetArrayF90(por,por_p,ierr)
    call VecGetArrayF90(rkin,rkin_p,ierr)
    call VecGetArrayF90(rrkin,rrkin_p,ierr)
  endif

  end subroutine ptran_init
  
  subroutine ptran_chem(da,da_1dof,da_kin,sles)

  use ptran_global_module
  use trdynmem_module
  use ptran_speciation_module
  use trgamdh_module
  use ptran_psi_module
  use ptran_setbnd_module

  implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#ifdef USE_PETSC216
    ! petsc-2.1.6
#include "include/finclude/petscsles.h"
#endif
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA    :: da, da_1dof, da_kin

#ifdef USE_PETSC216
  SLES  :: sles
#endif

#ifdef USE_PETSC221
  integer :: sles
#endif

! PetscViewer   vecviewer
  integer :: ierr,i,j,l,m
! integer :: k,ir,ii1,ii2,jj1,jj2,kk1,kk2,n,ng,nr
  integer :: isumjl
  real*8, pointer :: ccloc_p(:),temploc_p(:),ssat_loc_p(:),pressloc_p(:)

  if (myrank == 0) write(*,*) '--> initialize field variables'
  call trinit

  if (myrank == 0) write(*,*) '--> speciate initial fluid composition'
  call trstartup (da,da_1dof,da_kin,sles)

! initialize psi, gam, gamx
   
  call DAGlobalToLocalBegin(da_1dof,temp,INSERT_VALUES,temploc,ierr)
  call DAGlobalToLocalEnd(da_1dof,temp,INSERT_VALUES,temploc,ierr)

  call DAGlobalToLocalBegin(da,cc,INSERT_VALUES,ccloc,ierr)
  call DAGlobalToLocalEnd(da,cc,INSERT_VALUES,ccloc,ierr)

  call DAGlobalToLocalBegin(da_1dof,ssat,INSERT_VALUES,ssat_loc,ierr)
  call DAGlobalToLocalEnd(da_1dof,ssat,INSERT_VALUES,ssat_loc,ierr)

  call DAGlobalToLocalBegin(da_1dof,sat,INSERT_VALUES,sat_loc,ierr)
  call DAGlobalToLocalEnd(da_1dof,sat,INSERT_VALUES,sat_loc,ierr)
  
  call DAGlobalToLocalBegin(da_1dof,press,INSERT_VALUES,pressloc,ierr)
  call DAGlobalToLocalEnd(da_1dof,press,INSERT_VALUES,pressloc,ierr)

  call VecGetArrayF90(temploc,temploc_p,ierr)
  call VecGetArrayF90(ccloc,ccloc_p,ierr)
  call VecGetArrayF90(ssat_loc,ssat_loc_p,ierr)
  call VecGetArrayF90(pressloc,pressloc_p,ierr)
  
  gam = one
  gamx = one

!   print *,"Ptran_chem:", ssat_loc_p


  
  if (iact == 1) call trgamdh (ccloc_p,temploc_p)
  call trpsi (ccloc_p,pressloc_p,temploc_p,ssat_loc_p)
  call VecRestoreArrayF90(temploc,temploc_p,ierr)
  call VecRestoreArrayF90(ccloc,ccloc_p,ierr)
  call VecRestoreArrayF90(ssat_loc,ssat_loc_p,ierr)
  call VecRestoreArrayF90(pressloc,pressloc_p,ierr)
! initialize total concentrations
  psi = ppsi
  if (iphase == 2) psig = ppsig
 ! print *,"Ptran_chem:"
 ! print *, psig
 ! print *, ppsig
  

! initialize boundary conditions and source/sinks
  if (myrank==0) write(*,*) '--> speciate boundary/source fluid composition'
  call trsetbnd

  allocate(dfill(ncomp*ncomp))
  allocate(ofill(ncomp*ncomp))
  allocate(irow1(ncomp))
  allocate(irow2(ncomp))
  allocate(icol1(ncomp))
  allocate(icol2(ncomp))
! if (iblkfmt == 0) then
    do j = 1, ncomp
      do l = 1, ncomp
        dfill(j+(l-1)*ncomp) = 0
        ofill(j+(l-1)*ncomp) = 0
        isumjl = 0
        do i = 1, ncmplx
          isumjl = isumjl + shom(j,i)*shom(l,i)
        enddo
        if (isumjl.ne.0 .or. j.eq.l) then
          dfill(j+(l-1)*ncomp) = 1
          ofill(j+(l-1)*ncomp) = 1
        endif
        isumjl = 0
        do m = 1, nkin
          isumjl = isumjl + skin(j,m)*skin(l,m)
        enddo
        if (isumjl .ne. 0) then
          dfill(j+(l-1)*ncomp) = 1
        endif
      enddo
    enddo
    if (myrank == 0) then
      write(iunit2,'(/,"dfill:")')
      write(*,'(/,"dfill:")')
      do j = 1, ncomp
        write(iunit2,'(2x,30i2)') (dfill(j+(l-1)*ncomp),l=1,ncomp)
        write(*,'(2x,30i2)') (dfill(j+(l-1)*ncomp),l=1,ncomp)
      enddo
      write(iunit2,'(/,"ofill:")')
      write(*,'(/,"ofill:")')
      do j = 1, ncomp
        write(iunit2,'(2x,30i2)') (ofill(j+(l-1)*ncomp),l=1,ncomp)
        write(*,'(2x,30i2)') (ofill(j+(l-1)*ncomp),l=1,ncomp)
      enddo
    endif
    if (iblkfmt == 0) then
      call DASetBlockFills(da,dfill,ofill,ierr)
!     call DASetBlockFills(da,PETSC_NULL,ofill,ierr);CHKERRQ(ierr);
!     print *,'DASetBlockFills: ',ierr
    endif
! endif

  end subroutine ptran_chem

!===================================================================

  subroutine trinit
  
  use ptran_global_module
  use trdynmem_module
  use water_eos_module
  
  implicit none
  
  integer :: i,ii,j,iireg,jj,jaq,k,kk,kcount,kcountp,l,lp,m,mm,ml, &
             ns,ntemp,ncsq,ncsqdp,ncsqdpg,nr
  
  real*8 :: den,atm,btm,bdtm,frc,dadebye,dbdebye,dbext,bextend0,fjl

!-----set activity coefficient algorithm
      ntemp = ntmpmx
      i = 1
      do ii=1,ntemp
      
!       write(*,*) 'trinit: ',ii,tempini,temptab(ii)
        
        if(tempini.le.temptab(ii)) then
          i = ii
          goto 10
        endif
      enddo
      if (myrank == 0) write(iunit2,1001) tempini
 1001 format(' illegal temperature (',1pg12.4,') deg.c')
      if(tempini.le.300.) then
        if (myrank == 0) &
        write(*,*) 'error in temperature interpolation! STOP'
        stop
      endif

   10 continue

      if (i .eq. 1) then
        adebye  = atab(i)
        bdebye  = btab(i)
        bextend0 = bdotab(i)
      else
        atm = atab(i-1)
        btm = btab(i-1)
        bdtm = bdotab(i-1)
        frc   = (tempini-temptab(i-1))/(temptab(i)-temptab(i-1))
        dadebye = atab(i)   - atm
        dbdebye = btab(i)   - btm
        dbext   = bdotab(i) - bdtm
        adebye  = frc*dadebye + atm
        bdebye  = frc*dbdebye + btm
        bextend0 = frc*dbext  + bdtm
      endif
      if (iact.eq.1) then
        do j = 1, ncomp
          bextend(j) = bextend0
        enddo
        do i = 1, ncmplx
          bextendx(i) = bextend0
        enddo
      else if (iact .eq. 2) then ! Davies algorithm
        adebye = half
        bdebye = one
        do j = 1, ncomp
          a0(j) = one
          bextend(j) = 0.3d0*z(j)*z(j)*half
        enddo
        do i = 1, ncmplx
          ax0(i) = one
          bextendx(i) = 0.3d0*zx(i)*zx(i)*half
        enddo
      else if (iact .eq. 5) then
        iact = 1
        bextend0 = zero
        do j = 1, ncomp
          bextend(j) = bextend0
        enddo
        do i = 1, ncmplx
          bextendx(i) = bextend0
        enddo
      endif

!-----set initial temperature and density fields
      if (mode .eq. 2) then
!       do n = 1,nmax
!         temp(n) = tempini
!         if (pref0 .gt. zero) then
!           press(n) = pref0
!           call density(temp(n),pref0,den)
!           rho(n) = den*1.d-3
!         else
!           press(n) = -pref0
!           rho(n) = one
!         endif
!       enddo

        do ibc = 1,nblkbc
!         ibc = iregn(mmbc)
          if (pref0 .gt. zero) then
            call density (tempbc(ibc),pref0,den)
            dwbc(ibc) = den*1.d-3
          else
            dwbc(ibc) = 1.d0
          endif
        enddo
      endif
      
      ireg = initreg
      if (ibcreg(2).gt.0) ireg = ireg + ibcreg(2)-ibcreg(1) + 1
      if (ibcreg(4).gt.0) ireg = ireg + ibcreg(4)-ibcreg(3) + 1

      do iireg = 1, ireg
        do j = 1, ncomp
          if (itype(j,iireg) .eq. 0) then
            ncon(j,iireg) = 'log'
          else if (itype(j,iireg) .eq. 1) then
            ncon(j,iireg) = 'Total'
          else if (itype(j,iireg) .eq. 2) then
            ncon(j,iireg) = 'Aq+Sorp'
          else if (itype(j,iireg) .eq. 21) then
            ncon(j,iireg) = 'Aq+CEC'
          else if (itype(j,iireg) .eq. 5) then
            ncon(j,iireg) = 'Wt%'
          else if (itype(j,iireg) .eq. 7) then
            ncon(j,iireg) = 'Conc'
          else if (itype(j,iireg) .eq. 8) then
            ncon(j,iireg) = 'pH'
          else if (itype(j,iireg) .eq. 10) then
            ncon(j,iireg) = 'Constr. Qty'
          else if (itype(j,iireg) .eq. -1) then
            ncon(j,iireg) = 'Charge Bal.'
          endif
        enddo
      enddo

!-----set mineral kinetic rate law translation
      itypkini(20) = 1
      itypkini(21) = 2
      itypkini(22) = 3
      itypkini(25) = 4
      itypkini(30) = 5
      itypkini(90) = 6
      itypkini(91) = 7

!-----set aqueous kinetic rate law translation
      itypkiniaq(19) =  1
      itypkiniaq(20) =  2
      itypkiniaq(21) =  3
      itypkiniaq(22) =  4
      itypkiniaq(23) =  5
      itypkiniaq(25) =  6
      itypkiniaq(30) =  7
      itypkiniaq(60) =  8
      itypkiniaq(90) =  9
      itypkiniaq(91) = 10
      itypkiniaq(95) = 11
      
!-----load irreversible mineral properties
      do nr = 1, nkin
        do m = 1, mnrl
          if (namk(nr) .eq. namrl(m)) then
            ndxkin(nr) = m
            eqkin(nr) = alnk(m)
!           ze(nr) = ze0(m)
            do j = 1, ncpri
              skin(j,nr) = smnrl(j,m)
            enddo

!-----------store molar volume - liter/mole
            vbarkin(nr) = vbar(m)
            wtkin(nr)   = wtmin(m)
            goto 20
          endif
        enddo
        if (myrank == 0) write(*,*) 'mineral name not found: ',namk(nr)
        stop
   20   continue
   
!-----rate constant - mole/cm**2/sec

!-----convert units of kinetic reaction rate from mole/cm**3 to 
!     mole/liter
!     aeqkin=ten**eqkin(nr)

        if (rlim0(nr) .lt. zero) rlim0(nr) = ten**rlim0(nr)
        rlim(nr) = 1.d3*rlim0(nr)
        do l = npar1(nr), npar2(nr)
          if (rkf00(l) .lt. zero) rkf00(l) = ten**rkf00(l)
          rkf0(l)  = 1.d3*rkf00(l)
          rkf(l)   = rkf0(l)
        enddo
        if (rkfa00(nr) .lt. zero) rkfa00(nr) = ten**rkfa00(nr)
        if (rkfb00(nr) .lt. zero) rkfb00(nr) = ten**rkfb00(nr)
        rkfa0(nr) = rkfa00(nr)*1.d3
        rkfb0(nr) = rkfb00(nr)*1.d3
        rkfa(nr)  = rkfa0(nr)
        rkfb(nr)  = rkfb0(nr)
      enddo

!-----compress stoichiometric matrix storage
!-----homogeneous reactions

!-----store product of stoichiometric coefficients
      jj = 0
      do l = 1, ncomp
        do j = l, ncomp
          jj = jj + 1
          do i = 1, ncmplx
            sshom(i,jj) = shom(j,i)*shom(l,i)
          enddo
        enddo
      enddo

!-----compute nonzero dpsi(j,l) matrix indices: store in jpsi
      ncsqdp = 0
      jj = 0
      do l = 1, ncomp
        do j = 1, ncomp
          jj = jj + 1
          fjl = zero
          if (j.ne.l) then ! always include diagonal terms
            do i = 1, ncmplx
              fjl = shom(j,i)*shom(l,i)
              if (fjl .ne. zero) goto 555
            enddo
            goto 556
          endif
  555     continue
          ncsqdp = ncsqdp + 1
          jpsi(ncsqdp) = jj
  556     continue
        enddo
      enddo

!-----compute nonzero dpsig(j,l) matrix indices: store in jpsig
      if (ngas > 0 .and. iphase == 2) then
      ncsqdpg = 0
      jj = 0
      do l = 1, ncomp
        do j = 1, ncomp
          jj = jj + 1
          fjl = zero
          if (j.ne.l) then ! always include diagonal terms
            do i = 1, ngas
              fjl = sgas(j,i)*sgas(l,i)
              if (fjl .ne. zero) goto 557
            enddo
            goto 558
          endif
  557     continue
          ncsqdpg = ncsqdpg + 1
          jpsig(ncsqdpg) = jj
  558     continue
        enddo
      enddo
      endif
      
      ncsq = ncomp*ncomp
!     write(iunit2,*) 'trdatint-storage: ',ncsq,ncsqdp,ncsqdpg
!     do jj = 1, ncsqdp
!       jjj = jpsi(jj)
!       write(*,*) 'trdatint: ',jj,jjj,ncsqdp,ncomp*ncomp
!     enddo

      do j = 1, ncomp
        k = 0
        do i = 1, ncmplx
          if (shom(j,i).ne.zero) then
            k = k + 1
            fshom(j,k) = shom(j,i)
            ki(j,k) = i
          endif
        enddo
        lc(j) = k

        do l = j+1, ncomp
          kk = 0
          do i = 1, ncmplx
            if (shom(j,i).ne.zero .and. shom(l,i).ne.zero) then
              kk = kk + 1
              kki(j,l,kk) = i
              ffshom(j,l,kk) = shom(j,i)*shom(l,i)
            endif
          enddo
          llc(j,l) = kk
        enddo
      enddo

!-----gaseous reactions

!-----store product of stoichiometric coefficients
      if (ngas > 0 .and. iphase == 2) then
      k = 0
      do l = 1, ncomp
        do j = l, ncomp
          k = k + 1
          do i = 1, ngas
            ssgas(i,k) = sgas(j,i)*sgas(l,i)
          enddo
        enddo
      enddo

      do i = 1, ngas
        njg(i) = 0
        do j = 1, ncomp
          if (sgas(j,i).ne.zero) then
            njg(i) = njg(i) + 1
            jg(njg(i),i) = j
          endif
        enddo
      enddo

      do j = 1, ncomp
        k = 0
        do i = 1, ngas
          if (sgas(j,i).ne.zero) then
            k = k + 1
            fsgas(j,k) = sgas(j,i)
            kg(j,k) = i
          endif
        enddo
        lg(j) = k

        do l = j+1, ncomp
          kk = 0
          do i = 1, ngas
            if (sgas(j,i).ne.zero .and. sgas(l,i).ne.zero) then
              kk = kk + 1
              kkg(j,l,kk) = i
              ffsgas(j,l,kk) = sgas(j,i)*sgas(l,i)
            endif
          enddo
          llg(j,l) = kk
        enddo
      enddo
      endif
      
!-----mineral reactions
      do nr = 1, nkin
        k = 0
        do j = 1, ncomp
          if (skin(j,nr) .ne. zero) then
            k = k + 1
            kinj(nr,k) = j
            fskin(nr,k) = skin(j,nr)*skin(j,nr)
          endif
        enddo
        kmind(nr) = k
        k = 0
        do j = 1, ncomp
          do l = j+1, ncomp
            if (skin(j,nr).ne.zero .and. skin(l,nr).ne.zero) then
              k = k + 1
              kminj(nr,k) = j
              kminl(nr,k) = l
              ffskin(nr,k) = skin(j,nr)*skin(l,nr)
            endif
          enddo
        enddo
        kkmin(nr) = k
      enddo

      idebug = 0
      if (idebug > 0 .and. myrank == 0) then
        write(iunit2,*)
        write(iunit2,*) ' stoichiometrix matrix compression'
        kcount  = 0
        kcountp = 0
        do j = 1, ncomp
          write(iunit2,'(1x,a12,1x,"no. primary species: ",i4)') &
          nam(j),lc(j)
          write(iunit2,*) ' complex  stoichiometric coef.'
          do k = 1, lc(j)
            kcount = kcount + 1
            write(iunit2,'(a12,1x,1pe12.4)') namcx(ki(j,k)),fshom(j,k)
          enddo
          do l = j+1, ncomp
            write(iunit2,'("no. pairs: ",2x,2a12,1x,i4)') nam(j), &
            nam(l),llc(j,l)
            write(iunit2,*) ' complex  stoichiometric coef.'
            do k = 1, llc(j,l)
              kcountp = kcountp + 1
              write(iunit2,'(a12,1x,1pe12.4)') namcx(kki(j,l,k)), &
              ffshom(j,l,k)
            enddo
          enddo
        enddo

        write(iunit2,*) ' nonzero terms: ',kcount,' pairs: ',kcountp
        write(iunit2,*) ' nc x ncx = ',ncomp*ncmplx, &
        ', nc x nc x ncx = ',ncomp*ncomp*ncmplx
        write(iunit2,'()')
      endif

      do i = 1, ncmplx
        k = 0
        do j = 1, ncomp
          if (shom(j,i).ne.zero) then
            k = k + 1
            cshom(k,i) = shom(j,i)
            jcmpr(k,i) = j
          endif
        enddo
        ncmpr(i) = k
      enddo

!-----water stoichiometric coefficients
      if (myrank == 0) then
        write(iunit2,*)
        write(iunit2,*) 'stoichiometric coefficients for H2O'
        write(iunit2,*) 'species                   nH2O'
        if (jh2o.eq.0) then
          jaq = ncomp+1
        else
          jaq = jh2o
        endif
        do i = 1, ncmplx
          write(iunit2,'(1x,a20,1x,1pe12.4)') namcx(i),shom(jaq,i)
        enddo
        do i = 1, ngas
          write(iunit2,'(1x,a20,1x,1pe12.4)') namg(i),sgas(jaq,i)
        enddo
        do i = 1, nkin
          write(iunit2,'(1x,a20,1x,1pe12.4)') namk(i),skin(jaq,i)
        enddo
        do i = 1, nsrfmx
          write(iunit2,'(1x,a20,1x,1pe12.4)') namscx(i),ssorp(jaq,i)
        enddo
        write(iunit2,*)
      endif

!-----store mineral kinetic prefactor stoichiometry in proper order
      do m = 1, nkin
        do ml = npar1(m), npar2(m)

!---------set primary species prefactor coefficients
          do j = 1, nkinpri(ml)
            do l = 1, ncomp
              if (namprik(j,ml) .eq. nam(l)) then
                jpri(j,ml) = l
                goto 222
              endif
            enddo
            write(*,*) 'error finding prefactor primary species: stop', &
            myrank
            stop
  222       continue
          enddo

!---------set secondary species prefactor coefficients
          do i = 1, nkinsec(ml)
            do l = 1, ncmplx
              if (namseck(i,ml) .eq. namcx(l)) then
                isec(i,ml) = l
                goto 223
              endif
            enddo
            write(*,*) 'error finding prefactor secondary species: stop', &
            myrank
            stop
  223       continue
          enddo
        enddo
      enddo

   if (myrank == 0) then
      if (nkin > 0) then
        write(iunit2,'(/,"parallel reaction stoichiometry")')
        write(iunit2,'("mineral",10x,"      npri  nsec")')
      endif
      do nr = 1, nkin
        do lp = npar1(nr), npar2(nr)
          write(iunit2,'(a20,2i6)') namk(nr),nkinpri(lp),nkinsec(lp)
          do l = 1, nkinpri(lp)
            write(iunit2,'(10x,a20,1pe12.4)') nam(jpri(l,lp)), &
            skinpri(l,lp)
          enddo
          do l = 1, nkinsec(lp)
            write(iunit2,'(10x,a20,1pe12.4)') namcx(isec(l,lp)), &
            skinsec(l,lp)
          enddo
        enddo
      enddo
      write(iunit2,*)
   endif

!-----setup ion-exchange ordering of minerals, colloids, and cations
      ncollex = 0
      do m = 1, nexsolid
!-------set species indices for minerals
        do ns = 1, nkin
          if (namex(m) .eq. namk(ns)) then
            mex(m) = ns
            goto 50
          endif
        enddo

!-------set species indices for ion exchange colloids
        do j = ncomp-ncoll+1, ncomp
          if (namex(m) .eq. nam(j)) then
            mex(m) = j
            ncollex = ncollex + 1
            goto 50
          endif
        enddo

        write(*,*) 'error in trdatint setting ion exchange colloid ', &
        'and mineral indices: stop'
        stop

   50   continue

!-------set species indices for cations
        do jj = 1, nexmax
          do l = 1, ncomp
            if (namcat(jj) .eq. nam(l)) then
              jex(jj) = l
              goto 60
            endif
          enddo
          write(*,*) 'error in trdatint setting cation exchange ', &
          'indices: stop'
          stop
   60     continue
        enddo
      enddo

!-----determine nr. of colloid surface complexation and exchange solids
      ncolsrf = 0
      do m = 1, nsrfmin
        do mm = 1, ncoll
          if (namcoll(mm).eq.namsrf(m)) then
            ncolsrf = ncolsrf + 1
          endif
        enddo
      enddo
      
      do i = 1, nsrfmx
!       print *,i,nsrfmx,namscx(i),eqsorp0(i),eqsorp(i)
        if (eqsorp0(i).ne.zero) eqsorp(i) = eqsorp0(i)
      enddo

!-----setup mineral ordering for surface complexation reactions
      do mm = 1, nsrfmin
        do m = 1, nkin
          if (namsrf(mm) .eq. namk(m)) then
            msorp(mm) = m
            if (wtkin(m).le.zero .or. vbarkin(m).le.zero) then
              write(*,'(''error: zero formula weight or molar volume '', &
    &         ''for mineral '',a20,'': '',1p2e12.4)') namk(m),wtkin(m), &
              vbarkin(m)
              stop
            endif
            goto 70
          endif
        enddo

!-------set species indices for surface complexation colloids
        do j = ncomp-ncoll+1, ncomp
          if (namsrf(mm) .eq. nam(j)) then
            msorp(mm) = j
            if (wt(j).le.zero) then
              write(*,'(''error: zero formula weight '', &
    &         ''for colloid '',a20,'': '',1p2e12.4)') nam(j),wt(j)
              stop
            endif
            goto 70
          endif
        enddo

        write(*,*) 'error in trdatint: surface complex mineral not ', &
        'found! STOP',namsrf(mm),nam(j)
        stop
   70   continue
      enddo

  end subroutine trinit

end module ptran_init_module
