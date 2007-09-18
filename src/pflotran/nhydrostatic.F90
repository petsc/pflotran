
! This subroutine works only for stuctured grid 
module nhydrostatic_module

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

real*8, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

real*8, pointer :: hys_pres(:), hys_temp(:)

public nhydrostatic

contains

subroutine Get_Hydrosta_Pres(nz, dz, pref0, tref, dtdz, gravity, m_nacl) 
  use water_eos_module
  implicit none
  
  integer nz
  real*8 pref0, tref, dz(1:nz+1), dtdz, gravity, m_nacl
  real*8,save :: pref 
  data pref/-10.D0/ 
     
  integer n, itrho
  real*8 depth, pres0, tmp, rho0, rho1, rho, pres, dw_kg
  real*8 dzz,zz, xm_nacl, p 
! first go with z
      
   if(dabs(pref-pref0)<1D-10) return   
      
     xm_nacl = m_nacl * fmwnacl
     xm_nacl = xm_nacl /(1.D3 + xm_nacl)


      pref=pref0
      depth =0.D0
      pres0= pref
      tmp = tref
    
      
      rho0=rho
      p=pres0
      
    
    
      call nacl_den(tmp, p*1D-6, xm_nacl, dw_kg) 
      rho0 = dw_kg * 1D3

    
      
        do n = 1, nz+1
           if (n.eq.1) then
            dzz = 0.5d0*dz(1)
             zz = dzz
           elseif(n.eq.nz+1)then 
             dzz = 0.5d0*dz(nz)
             zz= zz+ dzz 
           else
             dzz = 0.5d0*(dz(n)+ dz(n-1))
             zz = zz + dzz
      endif
      tmp = tref +  dTdz*zz
    
  !   call wateos(tmp, p, rho, dw_mol, dwp, &
  !    dum, dum, dum, dum, grid%scale, ierr)
      call nacl_den(tmp, p*1D-6, xm_nacl, dw_kg) 
      rho = dw_kg * 1D3
 
                  
        itrho= 0
          
              !betap = rho * grid%gravity * grid%beta
              if(n==1)then
                pres = p + rho0 * gravity * dzz !+  betap * horiz
         !        call wateos(tmp, pres, rho, dw_mol, dwp, &
         !               dum, dum, dum, dum, grid%scale, ierr)
                  call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
                  rho = dw_kg * 1D3
                  

               else
               do
                pres = p + (rho0*dz(n-1) + rho* dz(n))/(dz(n)+dz(n-1))&
                     * gravity * dzz 
              !  call wateos(tmp, pres, rho1, dw_mol, dwp, &
              !  dum, dum, dum, dum, grid%scale, ierr)
                 call nacl_den(tmp, pres*1D-6, xm_nacl, dw_kg) 
                 rho1 = dw_kg * 1D3
                if (abs(rho1-rho) < 1.d-10) exit
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
     hys_pres(n)=p; hys_temp(n)= tmp
!    print *, 'nhydro: ', n, p, tmp

   enddo
end subroutine  Get_Hydrosta_Pres   
     
    


subroutine nhydrostatic(grid)

  use water_eos_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  PetscScalar, pointer :: xx_p(:) 
  
! integer :: itrho
  integer  ierr,i,j,k
  real*8 :: depth,horiz,dw_kg,p,dp
!           dum,dwp,rho,rho1,dw_mol,zz,dzz,tmp,pres,
            

! real*8 :: ibndtyp(grid%nblkbc),cbc(grid%nblkbc),sbc(grid%nblkbc)
! real*8 :: xxbc_rec(grid%ndof,grid%nblkbc),iphasebc_rec(grid%nblkbc)
  real*8 :: rho_ref, xm_nacl, pref

! real*8 :: betap,dz1,dz2,dx1,dx2, 
  integer na,nz, nc, m, nl, ibc, ng
  
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

 allocate(hys_pres(grid%nz+1))
 allocate(hys_temp(grid%nz+1))


  ! this 
  if (grid%iread_init /= 1) then
    call VecGetArrayF90(grid%xx, xx_p, ierr)

    !do nl=1, grid%nlmax
    xm_nacl = grid%m_nacl * fmwnacl
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)
    call nacl_den(grid%tref, grid%pref*1D-6, xm_nacl, dw_kg) 
    rho_ref = dw_kg * 1D3

    do i=grid%istart,grid%iend
      do j=grid%jstart,grid%jend
        do k=grid%kstart,grid%kend
            
          nl = 1 + (i-grid%istart) + (j-grid%jstart)*grid%nlx + (k-grid%kstart)*grid%nlxy  
          na=grid%nL2A(nl)+1 ! the natural ordering start from 0
        ! Get nx,ny,nz
            
          nz= floor(((real(na)-.5))/grid%nxy) + 1
          if(nz ==11) print*, i,j,k,nl, na

          ng=grid%nL2G(nl)
!geh          depth = grid%z(na)
!geh          horiz = grid%x(na)
          depth = grid%z(ng)
          horiz = grid%x(ng)
!geh          dp = rho_ref * grid%gravity * grid%beta *  (grid%x(grid%nmax) - grid%x(na))
          dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                     grid%x(ng))
          pref =  grid%pref + dp 
          call Get_Hydrosta_Pres(grid%nz,grid%dz0,pref, grid%tref, grid%dtdz, grid%gravity, grid%m_nacl) 
          xx_p(1+ (nl-1)*grid%ndof)=hys_pres(nz) 
          xx_p(2+ (nl-1)*grid%ndof)=hys_temp(nz)
        ! print *, i,j,k,nl,na,nz, hys_pres(nz), hys_temp(nz)
        enddo
      enddo  
    enddo
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  endif
  
  ! boundary conditions

  p = grid%pref
  
  
  depth = depth + 0.5d0*grid%dz0(grid%nz)
  horiz = horiz + 0.5d0*grid%dx0(grid%nx)
  
  dp = rho_ref * grid%gravity * grid%beta * horiz
  
  if (grid%myrank == 0) then
    write(*,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ",1pe11.4, &
 &          "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
    write(IUNIT2,'(" --> hydrostatic: length= ",1pe11.4,"[m], depth= ", &
           & 1pe11.4, "[m], dp= ",1pe11.4,"[Pa]")') horiz,depth,dp
  endif
  
  !save input grid%ibndtyp values
  !note-bcon must be in the order: 3D-left, right, top, bottom, front, back
  !                                1D-top, bottom

  call VecGetArrayF90(grid%xx, xx_p, ierr)
  do nc = 1, grid%nconnbc
    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    na=grid%nL2A(m) +1 
   !  nz= floor(((real(na)-.5))/grid%nxy) + 1
   !  nx= floor(((real(na-(nz-1)*grid%nlzy)/grid%ny 1
  
    ibc= grid%ibconn(nc)
    if(grid%iface(ibc)/=3 .and. grid%iface(ibc)/=4)then
      select case(grid%ibndtyp(ibc))
        case(4)
      
        case(3)
          na=grid%nL2A(m) +1 
          ng = grid%nL2G(m)
          nz= floor(((real(na)-.5))/grid%nxy) + 1
          if(grid%iface(ibc) ==1)then
!geh            dp = rho_ref * grid%gravity * grid%beta * grid%x(grid%nmax) 
            dp = rho_ref * grid%gravity * grid%beta * (grid%x_max-0.5d0*grid%x(grid%nx))
          elseif(grid%iface(ibc) ==2)then
            dp=0.D0
          else   
!geh            dp = rho_ref * grid%gravity * grid%beta * (grid%x(grid%nmax) - grid%x(na))
            dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                       grid%x(ng))
          endif
          pref =  grid%pref + dp 
          call Get_Hydrosta_Pres(grid%nz,grid%dz0,pref, grid%tref, grid%dtdz, grid%gravity, grid%m_nacl)  
          grid%xxbc(1,nc) = hys_pres(nz)
          grid%xxbc(2:grid%ndof,nc) = xx_p(2+ (m-1)*grid%ndof: m*grid%ndof)
        case(2)
          grid%xxbc(:,nc) = xx_p(1+ (m-1)*grid%ndof: m *grid%ndof)
        case(1) 
          na=grid%nL2A(m) +1
          ng = grid%nL2G(m)
          nz= floor(((real(na)-.5))/grid%nxy) + 1
          if(grid%iface(ibc) ==1)then
!geh            dp = rho_ref * grid%gravity * grid%beta * grid%x(grid%nmax) 
            dp = rho_ref * grid%gravity * grid%beta * (grid%x_max-0.5d0*grid%dx0(grid%nx))
          elseif(grid%iface(ibc) ==2)then
            dp=0.D0
          else   
!geh            dp = rho_ref * grid%gravity * grid%beta * (grid%x(grid%nmax) - grid%x(na))
            dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                       grid%x(ng))
          endif
          pref =  grid%pref + dp
          call Get_Hydrosta_Pres(grid%nz,grid%dz0,pref, grid%tref, grid%dtdz, grid%gravity, grid%m_nacl)  
          grid%xxbc(1,nc) = hys_pres(nz)
          grid%xxbc(2,nc) = hys_temp(nz)
        end select 
      elseif( grid%iface(ibc)==4 ) then
        select case(grid%ibndtyp(ibc))
        case(4)
      
        case(3)
          nz = grid%nz
          !nx= mod(mod(na,grid%nxy),grid%nx) + 1
          ng = grid%nL2G(m)
!geh          dp = rho_ref * grid%gravity * grid%beta * (grid%x(grid%nmax) - grid%x(na))
          dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                     grid%x(ng))
          pref =  grid%pref + dp 
          call Get_Hydrosta_Pres(grid%nz,grid%dz0,pref, grid%tref, grid%dtdz, grid%gravity, grid%m_nacl)  
          grid%xxbc(1,nc) = hys_pres(nz+1)
          grid%xxbc(2:grid%ndof,nc) = xx_p(2+ (m-1)*grid%ndof: m *grid%ndof)
        case(2)
          grid%xxbc(:,nc) = xx_p(1+ (m-1)*grid%ndof: m*grid%ndof)
        case(1) 
          nz = grid%nz
          ng = grid%nL2G(m)
!geh          dp = rho_ref * grid%gravity * grid%beta * (grid%x(grid%nmax) - grid%x(na))
          dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                     grid%x(na))
          pref =  grid%pref + dp 
          call Get_Hydrosta_Pres(grid%nz,grid%dz0,pref, grid%tref, grid%dtdz, grid%gravity, grid%m_nacl)  
          grid%xxbc(1,nc) = hys_pres(nz+1)
          grid%xxbc(2,nc) = hys_temp(nz+1)
        end select 
      elseif( grid%iface(ibc)== 3 ) then
        select case(grid%ibndtyp(ibc))
        case(4)
      
        case(3,1)
          na=grid%nL2A(m)+1
          ng = grid%nL2G(m)
!geh          dp = rho_ref * grid%gravity * grid%beta * (grid%x(grid%nmax) - grid%x(na))
          dp = rho_ref * grid%gravity * grid%beta * ((grid%x_max-0.5d0*grid%dx0(grid%nx)) - &
                                                     grid%x(ng))
          pref =  grid%pref + dp 
          grid%xxbc(1,nc) = pref
          grid%xxbc(2:grid%ndof,nc) = grid%tref
        case(2)
          grid%xxbc(:,nc) = xx_p(1+ (m-1)*grid%ndof: m*grid%ndof)
      end select
    endif
 !  if(nz ==11) stop
 !  print *, m,na,nz, nc, ibc, grid%iface(ibc),grid%ibndtyp(ibc),  grid%xxbc(:,nc), xx_p(1+ (m-1)*grid%ndof: 2+ (m-1)*grid%ndof)

  enddo
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
 
  deallocate(hys_pres)
  deallocate(hys_temp)

  end subroutine nhydrostatic
end module nhydrostatic_module

