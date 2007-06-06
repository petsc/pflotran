
! This subroutine works only for stuctured grid 

subroutine Get_Hydrosta_Pres(nz,ndz,dz,pref, tref, p,t,grid,ierr) 
  use water_eos_module
  implicit none
  
  real*8 depth, pref, tref, p, t, dz(1:ndz), dx(1:ndx)
  integer ndz, nz, ierr
  
  !note: ndy and dy list not used for now.
  
    call wateos(tref, pref, rho, dw_mol, dwp, &
                dum, dum, dum, dum, grid%scale, ierr)

! first go with z
      
      depth =0.D0
      pres0=pref
      tmp = grdi%tref
      rho0=rho
        
      do n =1, nz
         if(n==1)
          depth =  0.5* dz(n)
         else
           depth = 0.5d0*(dz(n-1) + dz(n))
           weight=0.5 
         endif   
          tmp = grid%dTdz * depth + tmp
          weight = dz(n-1)/(dz(n-1) + dz(n))
         itrho= 0
        do 
       ! betap = rho * grid%gravity * grid%beta
        pres = ((1.D0-weight)*rho+ weight * rho0) * grid%gravity * depth  + pres0! - betap * horiz
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
    pres0=pres; rho0=rho 
   enddo
! then go with x direction
 
    do  n=1,ndx
       
       
    enddo   
    
end subroutine  Get_Hydrosta_Pres   
     
    


subroutine mmhydrostatic(grid)

  use water_eos_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  PetscScalar, pointer :: xx_p(:) 
  
  integer i,ibc,ibc0,ierr,itrho,ireg,j,jm,k,m,n,nl,na !,nx,ny,nz
  real*8 :: betap,depth,dz1,dz2,horiz,dx1,dx2, &
            dum,dwp,rho,rho1,dw_mol,zz,dzz,tmp,pres,p,dp

  real*8 :: ibndtyp(grid%nblkbc),cbc(grid%nblkbc),sbc(grid%nblkbc)
  real*8 :: xxbc_rec(grid%ndof,grid%nblkbc),iphasebc_rec(grid%nblkbc)

  
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


  ! this 
  if (grid%iread_init /= 1) then
    call VecGetArrayF90(grid%xx, xx_p, ierr)

    do nl=1, grid%nlmax
      na=grid%nL2A(nl)+1 ! the natural ordering start from 0
  !   Get nx,ny,nz
      nz= int(na/grid%nxy) + 1
      depth = grid%z(na)
      horiz = grid%x(na)
      
        
      call Get_Hydrosta_Pres(nz,grid%nz,grid%dz0,grid%pref, grid%tref, pres,tmp,grid,ierr) 

      xx_p(1+ (nl-1)*grid%ndof)=pres 
      xx_p(2+ (nl-1)*grid%ndof)=tmp
        
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



  do nc = 1, grid%nconnbc
     m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
      
    if(ng<=0)then
      print *, "Wrong boundary node index... STOP!!!"
      stop
     end if

     ibc= grid%ibconn(nc)
     if(grid%iface(ibc)/=3 .and. grid%iface(ibc)/=4)then
       select case(grid%ibndtyp(ibc))
        case(3)
          grid%xxbc(:,nc) = xx_p(1+ (m-1)*grid%ndof: 2+ (m-1)*grid%ndof)
        case(2)
          grid%xxbc(:,nc) = xx_p(1+ (m-1)*grid%ndof: 2+ (m-1)*grid%ndof)
        case(1) 
          na=grid%nL2A(m)+1
          nz= int(na/grid%nxy) + 1
          call Get_Hydrosta_Pres(nz,grid%nz,grid%dz0,grid%xxbc0(1,ibc), grid%xxbc0(2,ibc), pres,tmp,grid,ierr)  
          grid%xxbc(1,nc) = pres
          grid%xxbc(2,nc) = tmp
       end select 
     ! elseif( grid%iface(ibc)==4 .and. grid%ibndtyp(ibc) == ) then 
      endif
           
  enddo
 
 end subroutine mmhydrostatic


