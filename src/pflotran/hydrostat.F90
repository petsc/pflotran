module hydrostat_module

  use pflow_gridtype_module

private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

#include "definitions.h"

  PetscReal, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

  public :: hydrostatic, mhydrostatic, owghydrostatic,recondition_bc
 
contains

subroutine hydrostatic (grid)
 
  use water_eos_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  
  PetscInt :: i,ibc,ibc0,ierr,itrho,ireg,j,jm,k,m,n !,nx,ny,nz
  PetscReal :: betap,depth,dz1,dz2,horiz,dx1,dx2, rho0,&
            dum,dwp,rho,rho1,dw_mol,zz,dzz,tmp,pres,p,dp

  PetscReal :: ibndtyp(grid%nblkbc),cbc(grid%nblkbc),sbc(grid%nblkbc)

  Vec :: temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, temp4_nat_vec
  PetscReal, pointer :: pres_p(:), temp_p(:), conc_p(:), xmol_p(:), sat_p(:)

!#include "definitions.h"

! set up hydrostatic initial and boundary conditions
! dp/dz = rho * g, dp/dx = -rho*g*beta, dT/dz = dTdz
  
! dTdz = geothermal gradient                    ! C/m
! beta = lateral hydrostatic gradient           ! m/m [dh/dx = -mu vx/(k rho g)]
! tref = reference temperature at top of domain ! C
! pref = reference pressure at top of domain    ! Pa

! depth = depth of domain top from surface   ! m
! length = width of system domain ! m

! this implementation is only approximate---need to take into account
! density and viscosity change with p and T: rho = rho(T,p)

!  vz = -k/mu (dp/dz - rho g z) = 0
!  vx = -k/mu dp/dx, h = p/(rho g)

  if (grid%iread_init == 0) then
  
#define GEH_SCALABLE
#ifdef GEH_SCALABLE
    if (grid%use_2ph /= PETSC_TRUE) then
      call VecGetArrayF90(grid%pressure,pres_p,ierr)
      call VecGetArrayF90(grid%temp,temp_p,ierr)
      call VecGetArrayF90(grid%conc,conc_p,ierr)
      call VecGetArrayF90(grid%sat,sat_p,ierr)
    else
      call VecGetArrayF90(grid%pressure,pres_p,ierr)
      call VecGetArrayF90(grid%temp,temp_p,ierr)
      call VecGetArrayF90(grid%xmol,xmol_p,ierr)
      call VecGetArrayF90(grid%sat,sat_p,ierr)
    endif
#else
!begin removed by geh  - not scalable memory-wise; see solution below
    if(grid%use_2ph /= PETSC_TRUE) then
      call DACreateNaturalVector(grid%da_1_dof,temp1_nat_vec,ierr) ! pressure
      call VecDuplicate(temp1_nat_vec,temp2_nat_vec,ierr)  ! temperature
      call VecDuplicate(temp1_nat_vec,temp3_nat_vec,ierr)  ! concentration
      call DACreateNaturalVector(grid%da_nphase_dof,temp4_nat_vec,ierr) ! saturation
       print *,'hydro: Create End'
    else
      call DACreateNaturalVector(grid%da_nphase_dof,temp1_nat_vec,ierr) ! pressure
      call VecDuplicate(temp1_nat_vec,temp4_nat_vec,ierr)  ! saturation
      call VecDuplicate(temp1_nat_vec,temp3_nat_vec,ierr)  ! concentration
      call DACreateNaturalVector(grid%da_1_dof,temp2_nat_vec,ierr) ! Temperature
    endif
     
    if (grid%myrank == 0) then
#endif
      ! set initial pressure and temperature fields
      dz1 = 0.d0
      dz2 = grid%dz0(1)
      depth = 0.d0
      do k = 1, grid%nz
        depth = depth + 0.5d0*(dz1+dz2)
        dz1 = dz2
        if (k < grid%nz) dz2 = grid%dz0(k+1)
        do j = 1, grid%ny
          dx1 = 0.d0
          dx2 = grid%dx0(1)
          horiz = 0.d0
          do i = 1, grid%nx
            horiz = horiz + 0.5d0*(dx1+dx2)
            dx1 = dx2
            if (i < grid%nx) dx2 = grid%dx0(i+1)
            tmp = grid%dTdz * depth + grid%tref
            betap = rho * grid%gravity * grid%beta
            pres = rho * grid%gravity * depth + grid%pref - betap * horiz
        
            call wateos(tmp, pres, rho, dw_mol, dwp, &
            dum, dum, dum, dum, grid%scale, ierr)
            
            itrho= 0
            do 
              betap = rho * grid%gravity * grid%beta
              pres = rho * grid%gravity * depth + grid%pref - betap * horiz
              call wateos(tmp, pres, rho1, dw_mol, dwp, &
              dum, dum, dum, dum, grid%scale, ierr)
              if (abs(rho-rho1) < 1.d-6) exit
              rho = rho1
              itrho = itrho + 1
              if (itrho > 100) then
                print *,' no convergence in hydrostat-stop',itrho,rho1,rho
                stop
              endif
            enddo
        
           
#ifdef GEH_SCALABLE
            if (i > grid%nxs .and. i <= grid%nxe .and. &
                j > grid%nys .and. j <= grid%nye .and. &
                k > grid%nzs .and. k <= grid%nze) then
              m = i-grid%nxs+(j-grid%nys-1)*grid%nlx+(k-grid%nzs-1)*grid%nlxy
              if (grid%use_2ph /= PETSC_TRUE) then
                jm = grid%jgas + m*grid%nphase - 1
                pres_p(m) = pres
                temp_p(m) = tmp
                conc_p(m) = grid%conc0
                sat_p(jm) = 1.d0
              else
                jm= m*grid%nphase
                pres_p(jm) = pres
                pres_p(jm+1) = pres
                temp_p(m) = tmp
                xmol_p(jm) = grid%conc0
                xmol_p(jm+1) = grid%conc0
                sat_p(jm) = 1.d0
                sat_p(jm+1) = 0.d0
              endif
            endif
#else
            m = i + (j-1)*grid%nx + (k-1)*grid%nxy - 1
!           print *,'hydrostatic: ',m+1,i,k,horiz,depth,tmp,pres,rho

            if (grid%use_2ph /= PETSC_TRUE) then
              jm = grid%jgas + m*grid%nphase - 1
              ! pressure
              ! print *,'hydro: set Val', grid%jgas, grid%nphase,m,jm
              call VecSetValue(temp1_nat_vec,m,pres,INSERT_VALUES,ierr)
              ! temperature
              call VecSetValue(temp2_nat_vec,m,tmp,INSERT_VALUES,ierr)
              ! conc
              call VecSetValue(temp3_nat_vec,m,grid%conc0,INSERT_VALUES,ierr)
              ! sat
              call VecSetValue(temp4_nat_vec,jm,1.d0,INSERT_VALUES,ierr)
            else    
              jm= m*grid%nphase
              ! pressure
              call VecSetValue(temp1_nat_vec,jm,pres,INSERT_VALUES,ierr)
              ! temperature
              call VecSetValue(temp1_nat_vec,jm+1,pres,INSERT_VALUES,ierr) 
              call VecSetValue(temp2_nat_vec,m,tmp,INSERT_VALUES,ierr)
              ! xl
              call VecSetValue(temp3_nat_vec,jm,grid%conc0,INSERT_VALUES,ierr)
              ! xg
              call VecSetValue(temp3_nat_vec,jm+1,grid%conc0,INSERT_VALUES,ierr)
              ! sl
              call VecSetValue(temp4_nat_vec,jm,1.d0,INSERT_VALUES,ierr)
              ! sg
              call VecSetValue(temp4_nat_vec,jm+1,0.d0,INSERT_VALUES,ierr) 
            endif 
#endif
          enddo
        enddo
      enddo
  
      !set co2 layer
#if 0
      if (grid%use_2ph /= PETSC_TRUE) then
        print *, 'ERROR: CO2 layer can only be set in mph mode'
        stop
      endif
      tmp = 50.d0
      pres = 2.d7
      do k = 1, 5
        do j = 1, 1
          do i = 1, 64
#ifdef GEH_SCALABLE
            if (i > grid%nxs .and. i <= grid%nxe .and. &
                j > grid%nys .and. j <= grid%nye .and. &
                k > grid%nzs .and. k <= grid%nze) then
              m = i-grid%nxs+(j-grid%nys-1)*grid%nlx+(k-grid%nzs-1)*grid%nlxy
              jm= m*grid%nphase
              pres_p(jm) = pres
              pres_p(jm+1) = pres
              temp_p(m) = tmp
              xmol_p(jm) = grid%conc0
              xmol_p(jm+1) = grid%conc0
              sat_p(jm) = 0.15d0
              sat_p(jm+1) = 0.85d0
              endif
            endif
#else
            m = i + (j-1)*grid%nx + (k-1)*grid%nxy - 1
            jm = m*grid%nphase
            ! pressure
            call VecSetValue(temp1_nat_vec,jm,pres,INSERT_VALUES,ierr)
            ! temperature
            call VecSetValue(temp1_nat_vec,jm+1,pres,INSERT_VALUES,ierr) 
            call VecSetValue(temp2_nat_vec,m,tmp,INSERT_VALUES,ierr)
            ! xl
            call VecSetValue(temp3_nat_vec,jm,grid%conc0,INSERT_VALUES,ierr)
            ! xg
            call VecSetValue(temp3_nat_vec,jm+1,grid%conc0,INSERT_VALUES,ierr)
            ! sl
            call VecSetValue(temp4_nat_vec,jm,0.15d0,INSERT_VALUES,ierr)
            ! sg
            call VecSetValue(temp4_nat_vec,jm+1,0.85d0,INSERT_VALUES,ierr) 
!            call VecSetValue(grid%pressure,m,pres,INSERT_VALUES,ierr)
!            call VecSetValue(grid%temp,m,tmp,INSERT_VALUES,ierr)
!            call VecSetValue(grid%conc,m,1.d0,INSERT_VALUES,ierr)
#endif
          enddo
        enddo
      enddo
#endif
  
      !set melt glass pressure, temperature, and tracer
#if 0
      if (grid%use_2ph == PETSC_TRUE) then
        print *, 'ERROR: Melt glass cannot be set in mph mode'
        stop
      endif
      tmp = 150.d0
!      pres = 1.d7
      pres = 8396027.61
      do k = 41, 45
        do j = 1, 1
          do i = 11, 31
#ifdef GEH_SCALABLE
            if (i > grid%nxs .and. i <= grid%nxe .and. &
                j > grid%nys .and. j <= grid%nye .and. &
                k > grid%nzs .and. k <= grid%nze) then
              m = i-grid%nxs+(j-grid%nys-1)*grid%nlx+(k-grid%nzs-1)*grid%nlxy
              pres_p(m) = pres
              temp_p(m) = tmp
              conc_p(m) = grid%conc0
              endif
            endif
#else
            m = i + (j-1)*grid%nx + (k-1)*grid%nxy - 1
            call VecSetValue(grid%pressure,m,pres,INSERT_VALUES,ierr)
            call VecSetValue(grid%temp,m,tmp,INSERT_VALUES,ierr)
            call VecSetValue(grid%conc,m,1.d0,INSERT_VALUES,ierr)
#endif
          enddo
        enddo
      enddo
#endif

#ifdef GEH_SCALABLE
    if (grid%use_2ph /= PETSC_TRUE) then
      call VecRestoreArrayF90(grid%pressure,pres_p,ierr)
      call VecRestoreArrayF90(grid%temp,temp_p,ierr)
      call VecRestoreArrayF90(grid%sat,sat_p,ierr)
      call VecRestoreArrayF90(grid%conc,conc_p,ierr)
    else
      call VecRestoreArrayF90(grid%pressure,pres_p,ierr)
      call VecRestoreArrayF90(grid%temp,temp_p,ierr)
      call VecRestoreArrayF90(grid%sat,sat_p,ierr)
      call VecRestoreArrayF90(grid%xmol,xmol_p,ierr)
    endif
#else
    endif
    
    call VecAssemblyBegin(temp1_nat_vec,ierr)
    call VecAssemblyEnd(temp1_nat_vec,ierr)
    call VecAssemblyBegin(temp2_nat_vec,ierr)
    call VecAssemblyEnd(temp2_nat_vec,ierr)
    call VecAssemblyBegin(temp3_nat_vec,ierr)
    call VecAssemblyEnd(temp3_nat_vec,ierr)
    call VecAssemblyBegin(temp4_nat_vec,ierr)
    call VecAssemblyEnd(temp4_nat_vec,ierr)
    if (grid%use_2ph /= PETSC_TRUE) then
      call DANaturalToGlobalBegin(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
                                  grid%pressure,ierr)
      call DANaturalToGlobalEnd(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
                                grid%pressure,ierr)
      call DANaturalToGlobalBegin(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                  grid%temp,ierr)
      call DANaturalToGlobalEnd(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                grid%temp,ierr)
      call DANaturalToGlobalBegin(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                                  grid%conc,ierr)
      call DANaturalToGlobalEnd(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                                grid%conc,ierr)
      call DANaturalToGlobalBegin(grid%da_nphase_dof,temp4_nat_vec, &
                                  INSERT_VALUES,grid%sat,ierr)
      call DANaturalToGlobalEnd(grid%da_nphase_dof,temp4_nat_vec, &
                                INSERT_VALUES,grid%sat,ierr)
    else
      call DANaturalToGlobalBegin(grid%da_nphase_dof,temp1_nat_vec, &
                                  INSERT_VALUES,grid%pressure,ierr)
      call DANaturalToGlobalEnd(grid%da_nphase_dof,temp1_nat_vec, &
                                INSERT_VALUES,grid%pressure,ierr)
      call DANaturalToGlobalBegin(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                  grid%temp,ierr)
      call DANaturalToGlobalEnd(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                grid%temp,ierr)
      call DANaturalToGlobalBegin(grid%da_nphase_dof,temp3_nat_vec, &
                                  INSERT_VALUES,grid%xmol,ierr)
      call DANaturalToGlobalEnd(grid%da_nphase_dof,temp3_nat_vec, &
                                INSERT_VALUES,grid%xmol,ierr)
      call DANaturalToGlobalBegin(grid%da_nphase_dof,temp4_nat_vec, &
                                  INSERT_VALUES,grid%sat,ierr)
      call DANaturalToGlobalEnd(grid%da_nphase_dof,temp4_nat_vec, &
                                INSERT_VALUES,grid%sat,ierr)
    endif
    call VecDestroy(temp1_nat_vec,ierr)
    call VecDestroy(temp2_nat_vec,ierr)
    call VecDestroy(temp3_nat_vec,ierr)
    call VecDestroy(temp4_nat_vec,ierr)
#endif
  endif

! boundary conditions

  p = grid%pref
  call wateos(grid%tref, p, rho, dw_mol, dwp, &
              dum, dum, dum, dum, grid%scale, ierr)
              
!geh  depth = grid%z(grid%nmax) + 0.5d0*grid%dz0(grid%nz)
!geh  horiz = grid%x(grid%nmax) + 0.5d0*grid%dx0(grid%nx)
  depth = grid%z_max
  horiz = grid%x_max

  
  dp = rho * grid%gravity * grid%beta * horiz
  
  if (grid%myrank == 0) then
    write(*,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ",1pe11.4, &
 &          "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
    write(IUNIT2,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ", &
                 &  1pe11.4, "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
  endif
  
  !save input grid%ibndtyp values
  !note-bcon must be in the order: 3D-left, right, top, bottom, front, back
  !                                1D-top, bottom
  do ibc = 1, grid%nblkbc
    if ((grid%iface(ibc) .ne. ibc) .and. (grid%nx > 1)) then
      print *,'error in hydrostatic: bcon out of order ',grid%iface(ibc),ibc
      stop
    endif
    ibndtyp(ibc) = grid%ibndtyp(ibc)
    cbc(ibc) = grid%concbc0(ibc)
    sbc(ibc) = grid%sgbc0(ibc)
  enddo
  
  if (grid%nblkbc > 6) then
    print *,'error in bcon/hydrostatic: ',grid%nblkbc,ibndtyp
    print *,'nblkbc must be 6 or less in bcon keyword'
    stop
  endif

! left and right faces

  ibc = 0
  ibc0 = 0
  
  if (grid%nx > 1) then
  
  ! left face
  
  ibc0 = ibc0 + 1
  
  p=grid%pref + dp
  call wateos(grid%tref, p, rho0, dw_mol, dwp, &
    dum, dum, dum, dum, grid%scale, ierr)
  
  do n = 1, grid%nz
    ibc = ibc + 1
    
!   if (ibc .gt. MAXBCREGIONS) then
!     print *,'pflowgrid_setup: too many boundary condition regions-stop!'
!     stop
!   endif
    
    grid%iregbc1(ibc) = n
    grid%iregbc2(ibc) = n
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 1  
    grid%k1bc(ibc) = n
    grid%k2bc(ibc) = n
    grid%j1bc(ibc) = 1
    grid%j2bc(ibc) = grid%ny
    grid%i1bc(ibc) = 1
    grid%i2bc(ibc) = 1
    
    if (n.eq.1) then
      dzz = 0.5d0*grid%dz0(1)
      zz = dzz
    else
      dzz = 0.5d0*(grid%dz0(n)+grid%dz0(n-1))
      zz = zz + dzz
    endif
    tmp = grid%tref + grid%dTdz*zz
    
    call wateos(tmp, p, rho, dw_mol, dwp, &
        dum, dum, dum, dum, grid%scale, ierr)

              
    itrho= 0
         
             ! betap = rho * grid%gravity * grid%beta
          if(n==1)then
                pres = p + rho0 * grid%gravity * dzz !+  betap * horiz
          !      tmp=tmp = grid%tref + grid%dTdz*zz
                call wateos(tmp, pres, rho, dw_mol, dwp, &
                     dum, dum, dum, dum, grid%scale, ierr)

                else
               do  
                pres = p + (rho0*grid%dz0(n-1) + rho*grid%dz0(n))/(grid%dz0(n)+grid%dz0(n-1))&
                     * grid%gravity * dzz 
                
              call wateos(tmp, pres, rho1, dw_mol, dwp, &
              dum, dum, dum, dum, grid%scale, ierr)
              if (abs(rho-rho1) < 1.d-6) exit
              rho = rho1
              itrho = itrho + 1
              if (itrho > 100) then
                print *,' no convergence in hydrostat-stop',itrho,rho1,rho
                stop
              endif
            enddo
           endif
     p = pres
     rho0=rho
 !   p = grid%pref + rho*grid%gravity*zz + dp


!   print *,'left: ',n,zz,tmp,p,rho,dp

    grid%pressurebc0(1,ibc) = p
    grid%tempbc0(ibc) = tmp
    grid%concbc0(ibc) = cbc(ibc0)
    grid%sgbc0(ibc) = sbc(ibc0)
    grid%velocitybc0(1,ibc) = 0.d0
    
!   print *,'BC: ',grid%myrank,n,ibc,zz,grid%z(grid%nmax),grid%x(grid%nmax), &
!   grid%tempbc(ibc),grid%tempbc(ibc+grid%nz), &
!   grid%pressurebc(1,ibc),grid%pressurebc(1,ibc+grid%nz)
  enddo

  ! right face
  
  ibc0 = ibc0 + 1

  p = grid%pref

  call wateos(grid%tref, p, rho0, dw_mol, dwp, &
    dum, dum, dum, dum, grid%scale, ierr)

  do n = 1, grid%nz
    ibc = ibc + 1
    
    grid%iregbc1(ibc) = grid%nz+n
    grid%iregbc2(ibc) = grid%nz+n
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 2
    grid%k1bc(ibc) = n
    grid%k2bc(ibc) = n
    grid%j1bc(ibc) = 1
    grid%j2bc(ibc) = grid%ny
    grid%i1bc(ibc) = grid%nx
    grid%i2bc(ibc) = grid%nx
    
    if (n.eq.1) then
      dzz = 0.5d0*grid%dz0(1)
      zz = dzz
    else
      dzz = 0.5d0*(grid%dz0(n)+grid%dz0(n-1))
      zz = zz + dzz
    endif
    tmp = grid%tref + grid%dTdz*zz

    call wateos(tmp, p, rho, dw_mol, dwp, &
    dum, dum, dum, dum, grid%scale, ierr)

!   call cowat (tmp, p, rw, uw, ierr)

!    p = grid%pref + rho*grid%gravity*zz !- dp
!   p = p0 + rho*grid%gravity*zz

!   print *,'right: ',n,zz,tmp,p,p0,rho,dp,grid%nx,grid%ny,grid%nz


      itrho= 0
                      ! betap = rho * grid%gravity * grid%beta
              if(n==1)then
                pres = p + rho0 * grid%gravity * dzz !+  betap * horiz
               call wateos(tmp, pres, rho, dw_mol, dwp, &
                     dum, dum, dum, dum, grid%scale, ierr)
 
              else
                 do 
  
                 pres = p + (rho0*grid%dz0(n-1) + rho* grid%dz0(n))/(grid%dz0(n)+grid%dz0(n-1))&
                     * grid%gravity * dzz 
              
             
               call wateos(tmp, pres, rho1, dw_mol, dwp, &
              dum, dum, dum, dum, grid%scale, ierr)
              if (abs(rho-rho1) < 1.d-6) exit
              rho = rho1
              itrho = itrho + 1
              if (itrho > 100) then
                print *,' no convergence in hydrostat-stop',itrho,rho1,rho
                stop
              endif
            enddo
           endif   
     p = pres
     rho0=rho

    grid%pressurebc0(1,ibc) = p
    grid%tempbc0(ibc) = tmp
    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
    grid%sgbc0(ibc) = sbc(ibc0)
    grid%velocitybc0(1,ibc) = 0.d0

!   print *,'BC: ',grid%myrank,n,ibc,zz,grid%z(grid%nmax),grid%x(grid%nmax), &
!   grid%tempbc(ibc),grid%tempbc(ibc+grid%nz), &
!   grid%pressurebc(1,ibc),grid%pressurebc(1,ibc+grid%nz)
  enddo
  
  endif
  
! top
  
  ibc0 = ibc0 + 1

  ibc = ibc + 1
  grid%iregbc1(ibc) = ibc
  grid%iregbc2(ibc) = ibc
! grid%ibndtyp(ibc) = 1 
  grid%ibndtyp(ibc) = ibndtyp(ibc0)
  grid%iface(ibc) = 3
  grid%k1bc(ibc) = 1
  grid%k2bc(ibc) = 1
  grid%j1bc(ibc) = 1
  grid%j2bc(ibc) = grid%ny
  grid%i1bc(ibc) = 1
  grid%i2bc(ibc) = grid%nx

  grid%pressurebc0(1,ibc) = grid%pref
  grid%tempbc0(ibc) = grid%tref
  grid%concbc0(ibc) = cbc(ibc0) !grid%conc0 0.85d0 !grid%conc0 !grid%concbc(3) !
  grid%sgbc0(ibc) = sbc(ibc0)
  grid%velocitybc0(1,ibc) = 0.d0
  
! bottom
  
  ibc0 = ibc0 + 1

  ibc = ibc + 1
  grid%iregbc1(ibc) = ibc
  grid%iregbc2(ibc) = ibc
! grid%ibndtyp(ibc) = 2 
  grid%ibndtyp(ibc) = ibndtyp(ibc0)
  grid%iface(ibc) = 4
  grid%k1bc(ibc) = grid%nz
  grid%k2bc(ibc) = grid%nz
  grid%j1bc(ibc) = 1
  grid%j2bc(ibc) = grid%ny
  grid%i1bc(ibc) = 1
  grid%i2bc(ibc) = grid%nx

! grid%pressurebc(1,ibc) = grid%pref
  grid%pressurebc0(1,ibc) = p !grid%pref + rho * grid%gravity * depth
  grid%tempbc0(ibc) = grid%tref + grid%dTdz * depth
  grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
  grid%sgbc0(ibc) = sbc(ibc0)
  grid%velocitybc0(1,ibc) = 0.d0

  if (grid%ny > 1) then
!     front
  
    ibc0 = ibc0 + 1

    ibc = ibc + 1
    grid%iregbc1(ibc) = ibc
    grid%iregbc2(ibc) = ibc
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 5
    grid%k1bc(ibc) = 1
    grid%k2bc(ibc) = grid%nz
    grid%j1bc(ibc) = 1
    grid%j2bc(ibc) = 1
    grid%i1bc(ibc) = 1
    grid%i2bc(ibc) = grid%nx

    grid%pressurebc0(1,ibc) = grid%pref
    grid%tempbc0(ibc) = grid%tref
    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
    grid%sgbc0(ibc) = sbc(ibc0)
    grid%velocitybc0(1,ibc) = 0.d0
  
!     back
  
    ibc0 = ibc0 + 1

    ibc = ibc + 1
    grid%iregbc1(ibc) = ibc
    grid%iregbc2(ibc) = ibc
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 6
    grid%k1bc(ibc) = 1
    grid%k2bc(ibc) = grid%nz
    grid%j1bc(ibc) = grid%ny
    grid%j2bc(ibc) = grid%ny
    grid%i1bc(ibc) = 1
    grid%i2bc(ibc) = grid%nx
    tmp = grid%tref + grid%dTdz * depth
    p = grid%pref + rho * grid%gravity * depth
    call wateos(tmp, p, rho, dw_mol, dwp, &
                dum, dum, dum, dum, grid%scale, ierr)

    grid%pressurebc0(1,ibc) = p
    grid%tempbc0(ibc) = tmp
    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
    grid%sgbc0(ibc) = sbc(ibc0)
    grid%velocitybc0(1,ibc) = 0.d0
  endif
  
  grid%nblkbc = ibc
    
!   do ibc = 1, grid%nblkbc
!     grid%pressurebc(2,ibc) = grid%pressurebc(1,ibc)
!     grid%velocitybc(2,ibc) = grid%velocitybc(1,ibc)
!   enddo

  if (grid%myrank == 0 ) then !.and. grid%iprint >= 3) then
    print *,'--> write out file pflow.bc'
    open(IUNIT3, file="pflow.bc", action="write", status="unknown")
    write(IUNIT3,'("BCON : nblkbc = ",i4)') grid%nblkbc
    do ibc = 1, grid%nblkbc
      write(IUNIT3,'(":ibndtyp  iface",/3x,i2,6x,i2)') grid%ibndtyp(ibc), &
      grid%iface(ibc)
      if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       p          T          sl          C")')
      else if (grid%ibndtyp(ibc) == 2) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       vl         T          sl          C")')
      endif
      do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
        if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    grid%pressurebc0(1,ireg),grid%tempbc0(ireg), &
                    1.d0-grid%sgbc0(ireg),grid%concbc0(ireg)
        else
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    grid%velocitybc0(1,ireg),grid%tempbc0(ireg), &
                    1.d0-grid%sgbc0(ireg),grid%concbc0(ireg)
        endif
      enddo
    enddo
    close(IUNIT3)
  endif

end subroutine hydrostatic
  
subroutine recondition_bc(grid)
  use water_eos_module

  implicit none
    type(pflowGrid), intent(inout) :: grid
     PetscReal :: p, pres, rho, rho1, rho0,  dum
     PetscInt :: nx,ny,nz, na,  nc, m, iln, ng, itrho, ierr,ibc  
     PetscReal  tmp, dzz,zz 
     PetscReal, pointer :: xx_p(:)
               
 ! only work for right boundary now                                 
     if (grid%ihydrostatic >= 1 .and. grid%nxe == grid%nx)then
       call VecGetArrayF90(grid%xx, xx_p, ierr)
       nx=grid%nx
        do ny=1, grid%ny
            p=grid%pref
             do nz=1, grid%nz   
               na= nx-1 + (ny-1)*grid%nx + (nz-1)*grid%nxy  
               do iln = 1, grid%nlmax
                 if(na == grid%nL2A(iln))then 
                    exit 
                 endif
               enddo    
               !tmp = xx_p((iln-1)*grid%ndof +2)
                       
               if(nz==1)then
                 p=grid%pref  
                 dzz=0.5d0*grid%dz0(1) 
                 zz=dzz
                 tmp = grid%tref
                  call wateos(tmp, p, rho0, dum, dum, &
                    dum, dum, dum, dum, grid%scale, ierr)
                  pres = p + rho0 * grid%gravity * dzz
                  rho=rho0  
               else  
                 dzz = 0.5d0*(grid%dz0(nz)+grid%dz0(nz-1))
                 zz = zz + dzz
                 
                 tmp = grid%tref + grid%dTdz*zz
                 itrho= 0
                
                do 
                   pres = p + (rho0*grid%dz0(nz-1) + rho* grid%dz0(nz))/(grid%dz0(nz)+grid%dz0(nz-1))&
                     * grid%gravity * dzz 
                                        
                   call wateos(tmp, pres, rho1, dum,dum, &
                      dum, dum, dum, dum, grid%scale, ierr)
                   if (abs(rho-rho1) < 1.d-6) exit
                   rho = rho1
                   itrho = itrho + 1
                   if (itrho > 100) then
                    print *,' no convergence in hydrostat-stop',itrho,rho1,rho
                    stop
                   endif
                 enddo
                endif
                
                p = pres
                rho0=rho
                
                do nc = 1, grid%nconnbc
                 ibc=  grid%ibconn(nc) 
                 if( grid%ibndtyp(ibc) == 3 .and. grid%iface(ibc)==2)then 
                    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
                    ng = grid%nL2G(m)
                    !print *, nx,ny, nz,na, p, tmp
                    if(na == grid%nL2A(m)) then
                         grid%pressurebc(1,nc) = p
                        ! print *, nc, p, tmp, nx,ny,nz,na,m
                         grid%tempbc(nc) = tmp
                         exit
                    endif
                  endif       
                 enddo
             enddo
            enddo  
        call VecRestoreArrayF90(grid%xx, xx_p, ierr)
        endif
 
            


end subroutine recondition_bc
  
 ! the coordinate in provided by pflow_compute_xyz, make sure this subroutine is called after that 
subroutine mhydrostatic(realization)

  use water_eos_module

  use Realization_module
  use Option_module
  use Grid_module
  use Structured_Grid_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Connection_module

  implicit none

  type(realization_type) :: realization
  PetscReal, pointer :: xx_p(:) 
  
  PetscInt :: ibc,ibc0,ierr,itrho,ireg,iz !,nl,ng
! PetscInt :: i,j,jm,k,m,
  PetscReal :: betap,depth,horiz,dx1,dx2,rho0,&
            rho,rho1,zz,dzz,tmp,pres,p,dp

!  PetscReal :: ibndtyp(grid%nblkbc) !,cbc(grid%nblkbc),sbc(grid%nblkbc)
!  PetscReal :: xxbc_rec(grid%ndof,grid%nblkbc),iphasebc_rec(grid%nblkbc)
  PetscReal :: xm_nacl, dw_kg
  PetscInt :: natural_id,nz, local_id, ghosted_id, iconn
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
! PetscReal :: dz1,dz2,dum,dwp,dw_mol,
  !Vec :: temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, temp4_nat_vec

!#include "definitions.h"

! set up hydrostatic initial and boundary conditions
! dp/dz = rho * g, dp/dx = -rho*g*beta, dT/dz = dTdz
  
! dTdz = geothermal gradient                    ! C/m
! beta = lateral hydrostatic gradient           ! m/m [dh/dx = -mu vx/(k rho g)]
! tref = reference temperature at top of domain ! C
! pref = reference pressure at top of domain    ! Pa

! depth = depth of domain top from surface   ! m
! length = width of system domain ! m

! this implementation is only approximate---need to take into account
! density and viscosity change with p and T: rho = rho(T,p)

!  vz = -k/mu (dp/dz - rho g z) = 0
!  vx = -k/mu dp/dx, h = p/(rho g)
   
  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.D3 + xm_nacl)
   ! call nacl_den(t, p*1D-6, xm_nacl, dw_kg) 
   !  dw_kg = dw_kg * 1D3
  

  if (option%iread_init /= 1) then
    call VecGetArrayF90(field%xx, xx_p, ierr)
! set initial pressure and temperature fields
    pres=option%pref     
    do local_id=1, grid%nlmax
!geh      na=grid%nL2A(nl)+1 ! the natural ordering start from 0
!geh      depth = grid%z(na)
!geh      horiz = grid%x(na)
      ghosted_id=grid%nL2G(local_id)
      depth = dabs(grid%z(ghosted_id))
      horiz = grid%x(ghosted_id)
      dx1 = dx2
      !print *,'mhydro', nl,na,depth,horiz
      tmp = option%dTdz * depth + option%tref
      betap = rho * option%gravity(3) * option%beta
    
     ! call wateos(tmp, pres, rho, dw_mol, dwp, &
     !             dum, dum, dum, dum, grid%scale, ierr)
      
      call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
      rho = dw_kg * 1D3
      
        pres = rho * option%gravity(3) * depth + option%pref - betap * horiz

      itrho= 0
      do 
        betap = rho * option%gravity(3) * option%beta
        pres = rho * option%gravity(3) * depth + option%pref - betap * horiz
        !call wateos(tmp, pres, rho1, dw_mol, dwp, &
        !            dum, dum, dum, dum, grid%scale, ierr)
        call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
         rho1 = dw_kg * 1D3
  
        if (abs(rho-rho1) < 1.d-10) exit
        rho = rho1
        itrho = itrho + 1
        if (itrho > 100) then
          print *,' no convergence in hydrostat-stop',itrho,rho1,rho
          stop
        endif
      enddo

      xx_p(1+ (local_id-1)*option%ndof)=pres 
      xx_p(2+ (local_id-1)*option%ndof)=tmp
        
    enddo
    call VecRestoreArrayF90(field%xx, xx_p, ierr)
  endif
  
  ! boundary conditions

  p = option%pref
!  call wateos(grid%tref, p, rho, dw_mol, dwp, &
!              dum, dum, dum, dum, grid%scale, ierr)
  call nacl_den(option%tref, p*1D-6, xm_nacl, dw_kg) 
  rho = dw_kg * 1D3
 
!geh  depth = grid%z(grid%nmax) + 0.5d0*grid%dz0(grid%nz)
!geh  horiz = grid%x(grid%nmax) + 0.5d0*grid%dx0(grid%nx)

  depth = grid%z_max
  horiz = grid%x_max
  
  dp = rho * option%gravity(3) * option%beta * horiz
  
  if (option%myrank == 0) then
    write(*,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ",1pe11.4, &
 &          "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
    write(IUNIT2,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ", &
           & 1pe11.4, "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
  endif
  
  !save input grid%ibndtyp values
  !note-bcon must be in the order: 3D-left, right, top, bottom, front, back
  !                                1D-top, bottom

  boundary_condition => realization%boundary_conditions%first
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nL2A(local_id)+1 

      p = option%pref + dp
!      call wateos(grid%tref, p, rho0, dw_mol, dwp, &
!      dum, dum, dum, dum, grid%scale, ierr)
      call nacl_den(option%tref, p*1D-6, xm_nacl, dw_kg) 
      rho0 = dw_kg * 1D3

      do iz = 1, (local_id-1)/grid%structured_grid%nlxy+1

        if (iz.eq.1) then
          dzz = 0.5d0*grid%structured_grid%dz0(1)
          zz = dzz
        else
          dzz = 0.5d0*(grid%structured_grid%dz0(iz)+grid%structured_grid%dz0(iz-1))
          zz = zz + dzz
        endif
        tmp = option%tref + option%dTdz*zz
      
    !   call wateos(tmp, p, rho, dw_mol, dwp, &
    !    dum, dum, dum, dum, grid%scale, ierr)
        call nacl_den(tmp, p*1D-6, xm_nacl, dw_kg) 
        rho = dw_kg * 1D3
   
                    
        itrho= 0
            
        !betap = rho * grid%gravity * grid%beta
        if (iz == 1) then
          pres = p + rho0 * option%gravity(3) * dzz !+  betap * hori
         ! call wateos(tmp, pres, rho, dw_mol, dwp, &
         !             dum, dum, dum, dum, grid%scale, ierr)
          call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
          rho = dw_kg * 1D3
        else
          do
            pres = p + (rho0*grid%structured_grid%dz0(iz-1) +  &
                        rho*grid%structured_grid%dz0(iz))/ &
                        (grid%structured_grid%dz0(iz)+grid%structured_grid%dz0(iz-1)) * &
                                                             option%gravity(3) * dzz 
           ! call wateos(tmp, pres, rho1, dw_mol, dwp, &
           ! dum, dum, dum, dum, grid%scale, ierr)
            call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
            rho1 = dw_kg * 1D3
            if (abs(rho-rho1) < 1.d-10) exit
            rho = rho1
            itrho = itrho + 1
            if (itrho > 100) then
              print *,' no convergence in hydrostat-stop',itrho,rho1,rho
              stop
            endif
          enddo
        endif
        p = pres
        rho0=rho
      enddo
      boundary_condition%aux_real_var(1,iconn) = p
      boundary_condition%aux_real_var(2,iconn) = tmp
    enddo
    boundary_condition => boundary_condition%next
  enddo
#if 0
  if (grid%myrank == 0 ) then !.and. grid%iprint >= 3) then
    print *,'--> write out file pflow.bc'
    open(IUNIT3, file="pflow.bc", action="write", status="unknown")
    write(IUNIT3,'("BCON : nblkbc = ",i4)') grid%nblkbc
    do ibc = 1, grid%nblkbc
      write(IUNIT3,'(":ibndtyp  iface",/3x,i2,6x,i2)') grid%ibndtyp(ibc), &
            grid%iface(ibc)
      if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       p          T   ", &
                     & "       sl          C")')
       else if (grid%ibndtyp(ibc) == 2) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       vl         T   ", &
                     & "       sl          C")')
      endif
      do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
        if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                        grid%i1bc(ireg),grid%i2bc(ireg), &
                        grid%j1bc(ireg),grid%j2bc(ireg), &
                        grid%k1bc(ireg),grid%k2bc(ireg), &
                        grid%xxbc0(:,ireg)
        else
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                        grid%i1bc(ireg),grid%i2bc(ireg), &
                        grid%j1bc(ireg),grid%j2bc(ireg), &
                        grid%k1bc(ireg),grid%k2bc(ireg), &
                        grid%velocitybc0(1,ireg),grid%xxbc0(:,ireg)
        endif
      enddo
    enddo
    close(IUNIT3)
  endif
#endif
end subroutine mhydrostatic



 ! the coordinate in provided by pflow_compute_xyz, make sure this subroutine is called after that 
subroutine owghydrostatic(grid)

  use water_eos_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  PetscReal, pointer :: xx_p(:) 
  
  PetscInt :: ibc,ibc0,ierr,itrho,ireg,n,nl,ng
! PetscInt :: i,j,jm,k,m
  PetscReal :: betap,depth,horiz,dx1,dx2, &
            dum,dwp,rho,rho1,dw_mol,zz,dzz,tmp,pres,p,dp

  PetscReal :: ibndtyp(grid%nblkbc) !,cbc(grid%nblkbc),sbc(grid%nblkbc)
  PetscReal :: xxbc_rec(grid%ndof,grid%nblkbc),iphasebc_rec(grid%nblkbc)

! PetscReal :: dz1,dz2
  !Vec :: temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, temp4_nat_vec

!#include "definitions.h"

! set up hydrostatic initial and boundary conditions
! dp/dz = rho * g, dp/dx = -rho*g*beta, dT/dz = dTdz
  
! dTdz = geothermal gradient                    ! C/m
! beta = lateral hydrostatic gradient           ! m/m [dh/dx = -mu vx/(k rho g)]
! tref = reference temperature at top of domain ! C
! pref = reference pressure at top of domain    ! Pa

! depth = depth of domain top from surface   ! m
! length = width of system domain ! m

! this implementation is only approximate---need to take into account
! density and viscosity change with p and T: rho = rho(T,p)

!  vz = -k/mu (dp/dz - rho g z) = 0
!  vx = -k/mu dp/dx, h = p/(rho g)

  if (grid%iread_init /= 1) then
    call VecGetArrayF90(grid%xx, xx_p, ierr)
! set initial pressure and temperature fields
    rho=1000.D0
    do nl=1, grid%nlmax
!geh      na=grid%nL2A(nl)+1 ! the natural ordering start from 0
!geh      depth = grid%z(na)
!geh      horiz = grid%x(na)
      ng=grid%nL2G(nl)
      depth = grid%z(ng)
      horiz = grid%x(ng)
      
      dx1 = dx2
      !print *,'mhydro', nl,na,depth,horiz
      tmp = grid%dTdz * depth + grid%tref
      betap = rho * grid%gravity * grid%beta
      pres = rho * grid%gravity * depth + grid%pref - betap * horiz

      call wateos(tmp, 2D7, rho, dw_mol, dwp, &
                  dum, dum, dum, dum, grid%scale, ierr)

      itrho= 0
      do 
        betap = rho * grid%gravity * grid%beta
        pres = rho * grid%gravity * depth + grid%pref - betap * horiz
        call wateos(tmp, 2D7, rho1, dw_mol, dwp, &
                    dum, dum, dum, dum, grid%scale, ierr)
        if (abs(rho-rho1) < 1.d-6) exit
        rho = rho1
        itrho = itrho + 1
        if (itrho > 100) then
          print *,' no convergence in hydrostat-stop',itrho,rho1,rho
          stop
        endif
      enddo

      xx_p(1+(nl-1)*grid%ndof) = pres 
      !xx_p(2+ (nl-1)*grid%ndof)=tmp
        
    enddo
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  endif
  
  ! boundary conditions

  p = grid%pref
  call wateos(grid%tref, p, rho, dw_mol, dwp, &
              dum, dum, dum, dum, grid%scale, ierr)
  
  depth = depth + 0.5d0*grid%dz0(grid%nz)
  horiz = horiz + 0.5d0*grid%dx0(grid%nx)
  
  dp = rho * grid%gravity * grid%beta * horiz
  
  if (grid%myrank == 0) then
    write(*,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ",1pe11.4, &
             &          "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
    write(IUNIT2,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ", &
             & 1pe11.4, "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
  endif
  
  !save input grid%ibndtyp values
  !note-bcon must be in the order: 3D-left, right, top, bottom, front, back
  !                                1D-top, bottom


  if (grid%nblkbc > 6) then
    print *,'error in bcon/hydrostatic: ',grid%nblkbc,ibndtyp
    print *,'nblkbc must be 6 or less in bcon keyword'
    stop
  endif

  xxbc_rec=grid%xxbc0
  iphasebc_rec=grid%iphasebc0
  
  do ibc = 1, grid%nblkbc
    if ((grid%iface(ibc) .ne. ibc) .and. (grid%nx > 1)) then
      print *,'error in hydrostatic: bcon out of order ',grid%iface(ibc),ibc
      stop
    endif
    ibndtyp(ibc) = grid%ibndtyp(ibc)
! grid%xxbc0(:,ibc)
!    cbc(ibc) = grid%concbc0(ibc)
!    sbc(ibc) = grid%sgbc0(ibc)
  enddo
  
 
! left and right faces

  ibc = 0
  ibc0 = 0
  
  if (grid%nx > 1) then
  
    ! left face
  
    ibc0 = ibc0 + 1
  
    do n = 1, grid%nz
      ibc = ibc + 1
    
!     if (ibc .gt. MAXBCREGIONS) then
!       print *,'pflowgrid_setup: too many boundary condition regions-stop!'
!       stop
!     endif
    
      grid%iregbc1(ibc) = n
      grid%iregbc2(ibc) = n
!     grid%ibndtyp(ibc) = 2 
      grid%ibndtyp(ibc) = ibndtyp(ibc0)
      grid%iface(ibc) = 1  
      grid%k1bc(ibc) = n
      grid%k2bc(ibc) = n
      grid%j1bc(ibc) = 1
      grid%j2bc(ibc) = grid%ny
      grid%i1bc(ibc) = 1
      grid%i2bc(ibc) = 1
      
      if (n.eq.1) then
        dzz = 0.5d0*grid%dz0(1)
        zz = dzz
      else
        dzz = 0.5d0*(grid%dz0(n)+grid%dz0(n-1))
        zz = zz + dzz
      endif
      tmp = grid%tref + grid%dTdz*zz
    
      call wateos(tmp, 2D7, rho, dw_mol, dwp, &
                  dum, dum, dum, dum, grid%scale, ierr)
    
      p = grid%pref + rho*grid%gravity*zz + dp
!     p = p0 + rho*grid%gravity*zz + dp

!     print *,'left: ',n,zz,tmp,p,rho,dp

!      grid%pressurebc0(1,ibc) = p
!      grid%tempbc0(ibc) = tmp
!      grid%concbc0(ibc) = cbc(ibc0)
!      grid%sgbc0(ibc) = sbc(ibc0)
      grid%velocitybc0(:,ibc) = 0.d0
      grid%xxbc0(1,ibc) = p
       !grid%xxbc0(2,ibc) = tmp
      grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
      grid%iphasebc0(ibc)=iphasebc_rec(ibc0)     
!     print *,'BC: ',grid%myrank,n,ibc,zz,grid%z(grid%nmax),grid%x(grid%nmax), &
!     grid%tempbc(ibc),grid%tempbc(ibc+grid%nz), &
!     grid%pressurebc(1,ibc),grid%pressurebc(1,ibc+grid%nz)
    enddo

    ! right face
  
    ibc0 = ibc0 + 1

    p = grid%pref
    do n = 1, grid%nz
      ibc = ibc + 1
    
      grid%iregbc1(ibc) = grid%nz+n
      grid%iregbc2(ibc) = grid%nz+n
!     grid%ibndtyp(ibc) = 2 
      grid%ibndtyp(ibc) = ibndtyp(ibc0)
      grid%iface(ibc) = 2
      grid%k1bc(ibc) = n
      grid%k2bc(ibc) = n
      grid%j1bc(ibc) = 1
      grid%j2bc(ibc) = grid%ny
      grid%i1bc(ibc) = grid%nx
      grid%i2bc(ibc) = grid%nx
    
      if (n.eq.1) then
        dzz = 0.5d0*grid%dz0(1)
        zz = dzz
      else
        dzz = 0.5d0*(grid%dz0(n)+grid%dz0(n-1))
        zz = zz + dzz
      endif
      tmp = grid%tref + grid%dTdz*zz

      call wateos(tmp, 2D7, rho, dw_mol, dwp, &
                  dum, dum, dum, dum, grid%scale, ierr)

!     call cowat (tmp, p, rw, uw, ierr)

      p = grid%pref + rho*grid%gravity*zz !- dp
!     p = p0 + rho*grid%gravity*zz

!     print *,'right: ',n,zz,tmp,p,p0,rho,dp,grid%nx,grid%ny,grid%nz

!      grid%pressurebc0(1,ibc) = p
!      grid%tempbc0(ibc) = tmp
!      grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
!      grid%sgbc0(ibc) = sbc(ibc0)
      grid%velocitybc0(:,ibc) = 0.d0
      grid%xxbc0(1,ibc) = p
      ! grid%xxbc0(2,ibc) = tmp
      grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
      grid%iphasebc0(ibc)=iphasebc_rec(ibc0)
!     print *,'BC: ',grid%myrank,n,ibc,zz,grid%z(grid%nmax),grid%x(grid%nmax), &
!     grid%tempbc(ibc),grid%tempbc(ibc+grid%nz), &
!     grid%pressurebc(1,ibc),grid%pressurebc(1,ibc+grid%nz)
    enddo
  
  endif
  
! top
  
  ibc0 = ibc0 + 1

  ibc = ibc + 1
  grid%iregbc1(ibc) = ibc
  grid%iregbc2(ibc) = ibc
! grid%ibndtyp(ibc) = 1 
  grid%ibndtyp(ibc) = ibndtyp(ibc0)
  grid%iface(ibc) = 3
  grid%k1bc(ibc) = 1
  grid%k2bc(ibc) = 1
  grid%j1bc(ibc) = 1
  grid%j2bc(ibc) = grid%ny
  grid%i1bc(ibc) = 1
  grid%i2bc(ibc) = grid%nx

!  grid%pressurebc0(1,ibc) = grid%pref
!  grid%tempbc0(ibc) = grid%tref
!  grid%concbc0(ibc) = cbc(ibc0) !grid%conc0 0.85d0 !grid%conc0 !grid%concbc(3) !
!  grid%sgbc0(ibc) = sbc(ibc0)
  grid%xxbc0(1,ibc) = grid%pref
  !  grid%xxbc0(2,ibc) = grid%tref
  grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
  grid%iphasebc0(ibc)=iphasebc_rec(ibc0)
  grid%velocitybc0(:,ibc) = 0.d0
! print *,'BC: top',grid%myrank,ibc,grid%xxbc0(:,ibc)
! bottom
  
  ibc0 = ibc0 + 1

  ibc = ibc + 1
  grid%iregbc1(ibc) = ibc
  grid%iregbc2(ibc) = ibc
! grid%ibndtyp(ibc) = 2 
  grid%ibndtyp(ibc) = ibndtyp(ibc0)
  grid%iface(ibc) = 4
  grid%k1bc(ibc) = grid%nz
  grid%k2bc(ibc) = grid%nz
  grid%j1bc(ibc) = 1
  grid%j2bc(ibc) = grid%ny
  grid%i1bc(ibc) = 1
  grid%i2bc(ibc) = grid%nx

!    grid%pressurebc(1,ibc) = grid%pref
!    grid%pressurebc0(1,ibc) = p !grid%pref + rho * grid%gravity * depth
!    grid%tempbc0(ibc) = grid%tref + grid%dTdz * depth
!    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
!    grid%sgbc0(ibc) = sbc(ibc0)
  grid%xxbc0(1,ibc) = p
!     grid%xxbc0(2,ibc) = tmp
  grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
  grid%iphasebc0(ibc)=iphasebc_rec(ibc0)
  grid%velocitybc0(:,ibc) = 0.d0
! print *,'BC: bot',grid%myrank,ibc,grid%xxbc0(:,ibc)

  if (grid%ny > 1) then
!     front
  
    ibc0 = ibc0 + 1

    ibc = ibc + 1
    grid%iregbc1(ibc) = ibc
    grid%iregbc2(ibc) = ibc
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 5
    grid%k1bc(ibc) = 1
    grid%k2bc(ibc) = grid%nz
    grid%j1bc(ibc) = 1
    grid%j2bc(ibc) = 1
    grid%i1bc(ibc) = 1
    grid%i2bc(ibc) = grid%nx

!    grid%pressurebc0(1,ibc) = grid%pref
!    grid%tempbc0(ibc) = grid%tref
!    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
!    grid%sgbc0(ibc) = sbc(ibc0)
    grid%xxbc0(1,ibc) = p
!     grid%xxbc0(2,ibc) = tmp
    grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
    grid%iphasebc0(ibc)=iphasebc_rec(ibc0) 
    grid%velocitybc0(:,ibc) = 0.d0
  
!     back
  
    ibc0 = ibc0 + 1

    ibc = ibc + 1
    grid%iregbc1(ibc) = ibc
    grid%iregbc2(ibc) = ibc
!   grid%ibndtyp(ibc) = 2 
    grid%ibndtyp(ibc) = ibndtyp(ibc0)
    grid%iface(ibc) = 6
    grid%k1bc(ibc) = 1
    grid%k2bc(ibc) = grid%nz
    grid%j1bc(ibc) = grid%ny
    grid%j2bc(ibc) = grid%ny
    grid%i1bc(ibc) = 1
    grid%i2bc(ibc) = grid%nx
    tmp = grid%tref + grid%dTdz * depth
    p = grid%pref + rho * grid%gravity * depth
    call wateos(tmp, 2D7, rho, dw_mol, dwp, &
                dum, dum, dum, dum, grid%scale, ierr)

!    grid%pressurebc0(1,ibc) = p
!    grid%tempbc0(ibc) = tmp
!    grid%concbc0(ibc) = cbc(ibc0) !grid%conc0
!    grid%sgbc0(ibc) = sbc(ibc0)
    grid%xxbc0(1,ibc) = p
   ! grid%xxbc0(2,ibc) = tmp
    grid%xxbc0(3,ibc) =xxbc_rec(3,ibc0)
    grid%iphasebc0(ibc)=iphasebc_rec(ibc0)
    grid%velocitybc0(:,ibc) = 0.d0
  endif
  
  grid%nblkbc = ibc
  
  if (ibc > MAXBCREGIONS) then
    print *,'error: increase MAXBCREGIONS: ', ibc,' > ',MAXBCREGIONS
    stop
  endif
    
! do ibc = 1, grid%nblkbc
!   grid%pressurebc(2,ibc) = grid%pressurebc(1,ibc)
!   grid%velocitybc(2,ibc) = grid%velocitybc(1,ibc)
! enddo

  if (grid%myrank == 0 ) then !.and. grid%iprint >= 3) then
    print *,'--> write out file pflow.bc'
    open(IUNIT3, file="pflow.bc", action="write", status="unknown")
    write(IUNIT3,'("BCON : nblkbc = ",i4)') grid%nblkbc
    do ibc = 1, grid%nblkbc
      write(IUNIT3,'(":ibndtyp  iface",/3x,i2,6x,i2)') grid%ibndtyp(ibc), &
      grid%iface(ibc)
      if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       p          T   ", &
                       & "       sl          C")')
      else if (grid%ibndtyp(ibc) == 2) then
        write(IUNIT3,'(": i1  i2  j1  j2  k1  k2       vl         T   ", &
                       & "       sl          C")')
      endif
      do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
        if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                        grid%i1bc(ireg),grid%i2bc(ireg), &
                        grid%j1bc(ireg),grid%j2bc(ireg), &
                        grid%k1bc(ireg),grid%k2bc(ireg), &
                        grid%xxbc0(:,ireg)
        else
          write(IUNIT3,'(6i4,1pe14.6,1p10e12.4)') &
                        grid%i1bc(ireg),grid%i2bc(ireg), &
                        grid%j1bc(ireg),grid%j2bc(ireg), &
                        grid%k1bc(ireg),grid%k2bc(ireg), &
                        grid%velocitybc0(1,ireg),grid%xxbc0(:,ireg)
        endif
      enddo
    enddo
    close(IUNIT3)
  endif

end subroutine owghydrostatic  
  
end module hydrostat_module
!======================================================================

